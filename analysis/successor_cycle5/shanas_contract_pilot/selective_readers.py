"""Source/metadata selection only. Never execute deposited author code.

Notebook outputs are skipped lexically without JSON-decoding their values.
CSV row fields outside an explicit whitelist are consumed without retention.
Headers are metadata and may be decoded in full. Numeric outcome cells are never
converted, printed, returned, or stored by the selective CSV iterator.
"""
import csv
import io
import json
import zipfile
import xml.etree.ElementTree as ET


class NotebookSourceOnly:
    def __init__(self, text):
        self.text, self.i = text, 0

    def ws(self):
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def string(self):
        self.ws()
        start = self.i
        if self.text[self.i] != '"':
            raise ValueError('Expected JSON string')
        self.i += 1
        while True:
            char = self.text[self.i]
            self.i += 1
            if char == '\\':
                self.i += 1
            elif char == '"':
                return json.loads(self.text[start:self.i])

    def skip(self):
        """Consume a JSON value structurally, without decoding its contents."""
        self.ws()
        if self.text[self.i] == '"':
            self.i += 1
            while True:
                char = self.text[self.i]
                self.i += 1
                if char == '\\':
                    self.i += 1
                elif char == '"':
                    return
        elif self.text[self.i] in '[{':
            opening = self.text[self.i]
            closing = ']' if opening == '[' else '}'
            self.i += 1
            self.ws()
            while self.text[self.i] != closing:
                if self.text[self.i] in ',:':
                    self.i += 1
                else:
                    self.skip()
                self.ws()
            self.i += 1
        else:
            while self.i < len(self.text) and self.text[self.i] not in ',]} \r\n\t':
                self.i += 1

    def cell(self):
        self.ws()
        if self.text[self.i] != '{':
            raise ValueError('Expected cell object')
        self.i += 1
        result = {}
        while True:
            self.ws()
            if self.text[self.i] == '}':
                self.i += 1
                return result
            key = self.string()
            self.ws()
            if self.text[self.i] != ':':
                raise ValueError('Expected colon')
            self.i += 1
            self.ws()
            start = self.i
            self.skip()
            if key in ['cell_type', 'source']:
                result[key] = json.loads(self.text[start:self.i])
            self.ws()
            if self.text[self.i] == ',':
                self.i += 1

    def read(self):
        self.ws()
        if self.text[self.i] != '{':
            raise ValueError('Expected notebook object')
        self.i += 1
        cells = None
        while True:
            self.ws()
            if self.text[self.i] == '}':
                break
            key = self.string()
            self.ws()
            if self.text[self.i] != ':':
                raise ValueError('Expected colon')
            self.i += 1
            self.ws()
            if key == 'cells':
                if self.text[self.i] != '[':
                    raise ValueError('Expected cells array')
                self.i += 1
                cells = []
                self.ws()
                while self.text[self.i] != ']':
                    cells.append(self.cell())
                    self.ws()
                    if self.text[self.i] == ',':
                        self.i += 1
                        self.ws()
                self.i += 1
            else:
                self.skip()
            self.ws()
            if self.text[self.i] == ',':
                self.i += 1
        if cells is None:
            raise ValueError('No cells array')
        return cells


def csv_header(stream, delimiter=','):
    return next(csv.reader(stream, delimiter=delimiter))


def selected_csv_rows(stream, allowed, delimiter=','):
    """RFC-style quoted field scanner; outside-whitelist values are not built."""
    header = csv_header(stream, delimiter)
    missing = set(allowed) - set(header)
    if missing:
        raise ValueError('Missing approved columns: ' + repr(sorted(missing)))
    chosen = {i: name for i, name in enumerate(header) if name in allowed}
    index, buf, values, quoted, at_start = 0, [], {}, False, True
    pending = None
    while True:
        char = pending if pending is not None else stream.read(1)
        pending = None
        if char == '':
            if quoted:
                raise ValueError('Unclosed CSV quote')
            if not at_start or index or values:
                if index in chosen:
                    values[chosen[index]] = ''.join(buf)
                if index + 1 != len(header):
                    raise ValueError('CSV column count mismatch')
                yield values
            return
        keep = index in chosen
        if quoted:
            if char == '"':
                next_char = stream.read(1)
                if next_char == '"':
                    if keep:
                        buf.append('"')
                else:
                    quoted = False
                    pending = next_char
            elif keep:
                buf.append(char)
        elif char == '"' and at_start:
            quoted = True
            at_start = False
        elif char == delimiter or char in '\r\n':
            if keep:
                values[chosen[index]] = ''.join(buf)
            buf, at_start = [], True
            if char == delimiter:
                index += 1
            else:
                if char == '\r':
                    next_char = stream.read(1)
                    if next_char != '\n':
                        pending = next_char
                if index + 1 != len(header):
                    raise ValueError('CSV column count mismatch')
                yield values
                index, values = 0, {}
        else:
            at_start = False
            if keep:
                buf.append(char)


NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def xlsx_shared_strings(zf, wanted):
    """Retain only shared string entries needed by selected cells."""
    result = {}
    if not wanted:
        return result
    with zf.open('xl/sharedStrings.xml') as stream:
        index = 0
        for event, element in ET.iterparse(stream, events=['end']):
            if element.tag.endswith('}si'):
                if index in wanted:
                    result[index] = ''.join(t.text or '' for t in element.findall('.//s:t', NS))
                element.clear()
                index += 1
    return result


def xlsx_header_schema(data):
    """Read sheet labels and each first row only; no data-row cell decoding."""
    result = []
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        workbook = ET.fromstring(zf.read('xl/workbook.xml'))
        relationships = ET.fromstring(zf.read('xl/_rels/workbook.xml.rels'))
        targets = {x.attrib['Id']: x.attrib['Target'] for x in relationships}
        for sheet in workbook.findall('.//s:sheet', NS):
            rel = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target = targets[rel]
            path = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            cells = []
            with zf.open(path) as stream:
                for event, element in ET.iterparse(stream, events=['end']):
                    if element.tag.endswith('}row'):
                        for cell in element.findall('s:c', NS):
                            v = cell.find('s:v', NS)
                            inline = cell.find('s:is', NS)
                            cells.append({'ref': cell.attrib['r'], 'type': cell.get('t'),
                                          'raw': v.text if v is not None else '',
                                          'inline': ''.join(inline.itertext()) if inline is not None else ''})
                        break
            requested = {int(c['raw']) for c in cells if c['type'] == 's'}
            strings = xlsx_shared_strings(zf, requested)
            headers = {c['ref']: strings[int(c['raw'])] if c['type'] == 's' else
                       c['inline'] if c['type'] == 'inlineStr' else c['raw'] for c in cells}
            result.append({'sheet_name': sheet.attrib['name'], 'xml_path': path, 'first_row': headers})
    return result


def xlsx_selected_metadata(data, allowed, sheet_name):
    schema = next(x for x in xlsx_header_schema(data) if x['sheet_name'] == sheet_name)
    by_column = {''.join(c for c in ref if c.isalpha()): name
                 for ref, name in schema['first_row'].items() if name in allowed}
    if set(by_column.values()) != set(allowed):
        raise ValueError('Approved metadata columns missing from workbook')
    rows, requested = [], set()
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        with zf.open(schema['xml_path']) as stream:
            for event, element in ET.iterparse(stream, events=['end']):
                if not element.tag.endswith('}row'):
                    continue
                if element.get('r') == '1':
                    element.clear()
                    continue
                chosen = {}
                for cell in element.findall('s:c', NS):
                    col = ''.join(c for c in cell.attrib['r'] if c.isalpha())
                    if col not in by_column:
                        continue
                    v = cell.find('s:v', NS)
                    inline = cell.find('s:is', NS)
                    if cell.get('t') == 's':
                        val = int(v.text)
                        requested.add(val)
                        chosen[by_column[col]] = ('shared', val)
                    elif cell.get('t') == 'inlineStr':
                        chosen[by_column[col]] = ('text', ''.join(inline.itertext()))
                    else:
                        chosen[by_column[col]] = ('text', v.text if v is not None else '')
                rows.append(chosen)
                element.clear()
        strings = xlsx_shared_strings(zf, requested)
    return [{name: strings[item[1]] if item[0] == 'shared' else item[1]
             for name, item in row.items()} for row in rows]


def self_test():
    fixture = '{"cells":[{"cell_type":"code","source":["a=1\\n"],"outputs":[{"data":{"value":123456.789}}]}],"metadata":{"skip":[1,2]}}'
    if NotebookSourceOnly(fixture).read() != [{'cell_type': 'code', 'source': ['a=1\n']}]:
        raise RuntimeError('Notebook selector failed')
    fixture = 'ID,sequence,protected\na,"ACG",123.4\r\nb,"A""C",456.7\n'
    expected = [{'ID': 'a', 'sequence': 'ACG'}, {'ID': 'b', 'sequence': 'A"C'}]
    if list(selected_csv_rows(io.StringIO(fixture), {'ID', 'sequence'})) != expected:
        raise RuntimeError('CSV selector failed')
    workbook = io.BytesIO()
    with zipfile.ZipFile(workbook, 'w') as zf:
        zf.writestr('xl/workbook.xml', '<workbook xmlns="'+NS['s']+'" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Metadata" sheetId="1" r:id="rId1"/></sheets></workbook>')
        zf.writestr('xl/_rels/workbook.xml.rels', '<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        zf.writestr('xl/worksheets/sheet1.xml', '<worksheet xmlns="'+NS['s']+'"><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>ID</t></is></c><c r="B1" t="inlineStr"><is><t>sequence</t></is></c><c r="C1" t="inlineStr"><is><t>protected</t></is></c></row><row r="2"><c r="A2"><v>7</v></c><c r="B2" t="inlineStr"><is><t>ACGT</t></is></c><c r="C2"><v>DO_NOT_DECODE_AS_NUMBER</v></c></row></sheetData></worksheet>')
    if xlsx_selected_metadata(workbook.getvalue(), {'ID', 'sequence'}, 'Metadata') != [{'ID':'7', 'sequence':'ACGT'}]:
        raise RuntimeError('XLSX selector failed')
    return 'synthetic selectors passed; no real outcomes'


if __name__ == '__main__':
    print(self_test())
