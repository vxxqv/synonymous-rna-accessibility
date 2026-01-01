#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Standard sense code. Enumeration is limited before recursion; no sampling.
const std::string bases="UCAG";
const std::string amino="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
bool paired(char a,char b) {
  return (a=='A'&&b=='U')||(a=='U'&&(b=='A'||b=='G'))||
         (a=='G'&&(b=='C'||b=='U'))||(a=='C'&&b=='G');
}
int main(int argc,char**argv) {
 try {
  if(argc<2) throw std::runtime_error("usage: exact_kernel PEPTIDE [i:j ...] (zero based)");
  std::map<char,std::vector<std::string>> code;
  int q=0;
  for(char a:bases)for(char b:bases)for(char c:bases) {
    std::string s; s+=a;s+=b;s+=c;code[amino[q++]].push_back(s);
  }
  std::string peptide=argv[1], sequence(3*peptide.size(),'N');
  std::vector<std::pair<int,int>> edges;
  std::vector<int> used(sequence.size(),0);
  for(int i=2;i<argc;i++) {
    std::string token=argv[i];auto colon=token.find(':');
    if(colon==std::string::npos)throw std::runtime_error("invalid pair");
    int a=std::stoi(token.substr(0,colon)),b=std::stoi(token.substr(colon+1));
    if(a>b)std::swap(a,b);
    if(a<0||b>=int(sequence.size())||a==b||used[a]++||used[b]++)throw std::runtime_error("not a matching");
    edges.emplace_back(a,b);
  }
  std::sort(edges.begin(),edges.end());
  std::uint64_t total=1;
  for(char a:peptide) {
    if(a=='*'||!code.count(a))throw std::runtime_error("invalid sense peptide");
    if(total>50000/code[a].size())throw std::runtime_error("Cartesian cap exceeded");
    total*=code[a].size();
  }
  std::function<void(std::size_t)> visit=[&](std::size_t k) {
    if(k==peptide.size()) {
      std::cout<<sequence<<'\t';
      for(auto [i,j]:edges)std::cout<<((sequence[i]=='A'||sequence[i]=='G')?'R':'Y');
      for(std::size_t i=0;i<peptide.size();i++)if(peptide[i]=='S')std::cout<<(sequence[3*i]=='A'?'a':'u');
      std::cout<<'\n';return;
    }
    for(const auto& codon:code[peptide[k]]) {
      sequence.replace(3*k,3,codon);bool ok=true;
      for(auto [i,j]:edges)if(j<int(3*k+3)&&!paired(sequence[i],sequence[j])){ok=false;break;}
      if(ok)visit(k+1);
    }
  };
  visit(0);
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
