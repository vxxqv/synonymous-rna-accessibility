"""Generate exact figure-source arrays for the general-graph theory pilot."""
import itertools
import json
from pathlib import Path
from graph_compiler import CODONS,valid,normalize,endpoint_certificate,weighted_projection

def complete_graph(peptide,edges):
    vertices=[''.join(z)for z in itertools.product(*(CODONS[a]for a in peptide))if valid(peptide,edges,''.join(z))]
    edges_index=[(i,j)for i,x in enumerate(vertices)for j,y in enumerate(vertices)if i<j and sum(a!=b for a,b in zip(x,y))==1]
    return vertices,edges_index

def main():
    output=Path(__file__).resolve().parents[2]/'results'/'expansion';output.mkdir(parents=True,exist_ok=True)
    peptide='LKK';edges=[(0,5),(0,8)];x='CUCAAGAAG';y='UUAAAAAAA'
    certificate=endpoint_certificate(peptide,edges,x,y);vs,es=complete_graph(peptide,edges)
    paths=[];orders=[]
    for order in itertools.permutations(certificate['changed_positions0']):
        z=list(x);path=[x]
        for p in order:
            z[p]=y[p]
            if not valid(peptide,edges,''.join(z)):break
            path.append(''.join(z))
        else:paths.append(path);orders.append(order)
    assert len(paths)==2
    branching={'peptide':peptide,'pair_edges_positions0':edges,'endpoints':[x,y],'certificate':certificate,'vertices':vs,'edges_vertex_indices':es,'shortest_paths':paths,'shortest_mutation_orders_positions0':orders,'dag_event_labels':['Leu first C->U','Leu wobble C->A','Lys1 wobble G->A','Lys2 wobble G->A'],'dag_positions0':[[1,0],[0,0],[2,-0.6],[2,0.6]]}
    peptide='LK';edges=[(0,5)];vs,es=complete_graph(peptide,edges);projection=weighted_projection(peptide,edges,'CUCAAG')
    groups={}
    for x in vs:groups.setdefault((x[0],x[5]),[]).append(x)
    projected=[{'paired_bases':list(k),'multiplicity':len(xs),'sequences':xs,'danger_bits':[int(k[0]=='C'),int(k[1]=='A')]}for k,xs in sorted(groups.items())]
    assert sorted(g['multiplicity']for g in projected)==[2,2,4]
    weighted={'peptide':peptide,'pair_edges_positions0':edges,'vertices':vs,'edges_vertex_indices':es,'projection':projection,'projected_states':projected}
    odds=[{'partners':r,'volume':str(4+2**(r+1)),'C_first_probability_numerator':1,'C_first_probability_denominator':1+2**(r-1),'shortest_paths':str(__import__('math').factorial(r)),'changed_events':r+2}for r in range(1,31)]
    payload={'branching_LKK':branching,'weighted_Leu_Lys':weighted,'star_family':odds,'interpretation':'Compatibility counts under uniform synonymous genotypes; no energetic, fitness or measured transition-rate meaning.'}
    (output/'general_graph_figure_sources.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'branching_vertices':len(branching['vertices']),'branching_edges':len(branching['edges_vertex_indices']),'branching_shortest_paths':len(paths),'weighted_vertices':len(vs),'projected_states':len(projected),'projected_multiplicities':sorted(g['multiplicity']for g in projected)}))

if __name__=='__main__':main()
