#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <mpi.h>  
#include <random>   

#define alfah 0.01
#define alfaj 0.01
#define tau 1000
#define exponent 1
#define lenMean 50000
#define lenWait 10000
#define lenObs 5000
using namespace std;

unsigned nt2index(char nucleotide){
  if(nucleotide=='A') return 0;
if(nucleotide=='a') return 0;
  if(nucleotide=='U') return 1;
  if(nucleotide=='u') return 1;
  if(nucleotide=='T') return 1;
  if(nucleotide=='C') return 2;
  if(nucleotide=='c') return 2;
  if(nucleotide=='G') return 3;
  if(nucleotide=='g') return 3;
  if(nucleotide=='-') return 4;
}

char index2nt(unsigned index){
  if(index==0) return 'A';
  if(index==1) return 'U';
  if(index==2) return 'C';
  if(index==3) return 'G';
  if(index==4) return '-';
}

class histo{
  public:
  void set(unsigned i, unsigned j,double x);
  double count[5][5];
  double mean[5][5];
  double sing[5];
  double sum[5];
  histo(){
    for(unsigned i=0;i<5;++i){ sing[i]=0; sum[i]=0; for(unsigned j=0;j<5;++j){mean[i][j]=0.0; count[i][j]=0.0;}}
  }
};
int main(int argc,char*argv[]){
	MPI_Init(&argc,&argv);
	int rank,np;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
   	std::ifstream input(argv[1]);
	std::string line;
	std::vector<string> seq;
	std::vector<string> names;
	int iseq=0, seqsize, len, seqid;
	while(getline(input, line)){
		seq.resize(iseq+1);
		seq[iseq]+=line;
		iseq++;
	}
	seqsize=seq.size();
	len=seq[0].length();
	for(unsigned i=0;i<seqsize;++i)
		assert(seq[i].length()==len);
	double identity[seqsize];
	for(unsigned is=0;is<seqsize;++is){
		identity[is]=1.0;
		for(unsigned it=0;it<seqsize;++it){
			seqid=0;
			for(unsigned in=0;in<len;++in){
				char nt1=seq[is][in];
				char nt2=seq[it][in];
				if(nt2index(nt1)==nt2index(nt2)) seqid++;
			}
			if((double)seqid/len>=0.9) identity[is]++; 
		}
	}
	double freq[len][5], Jcum, Jcumold, deltaH, t, L,lh; 
	unsigned fh[len], iold, inew,jold;
	L=(double)atof(argv[2]);
	lh=L;
	std::vector<histo> histogram(len*len);
	std::vector<histo> J(len*len);
	std::vector<histo> h(len);
	std::vector<histo> new_J(len*len);
	std::vector<histo> new_h(len);
	std::vector<histo> F(len*len);
	std::vector<histo> new_F(len*len);
	std::vector<histo> f(len);
	std::vector<histo> new_f(len);

	std::random_device rd {}; 
	std::mt19937 g{rd()}; 
    std::uniform_int_distribution<int> unif(0,seqsize-1);
	std::uniform_int_distribution<int> nucl(0,4);
	std::uniform_real_distribution<double> rando(0.0,1.0);

	//INIZIALIZZAZIONI	
	unsigned init=unif(g);
	for(unsigned in=0;in<len;in++){
		char nt1=seq[init][in];
		fh[in]=nt2index(nt1);
		for(unsigned i=0; i<5; i++) freq[in][i]=0;
	}	
	//ISTOGRAMMI ALLINEAMENTO
	for(unsigned in=0;in<len;in++){ for(unsigned jn=0;jn<len;jn++){
		for(unsigned is=0;is<seqsize;++is){
			char nt1=seq[is][in];
			char nt2=seq[is][jn];
			histogram[in*len+jn].count[nt2index(nt1)][nt2index(nt2)]+=(1.0/identity[is]);
			}
		}
		for(unsigned is=0;is<seqsize;++is){
			char nt1=seq[is][in];
			freq[in][nt2index(nt1)]+=(1.0/identity[is]);
	}}
	double Meff=0.0;
        double kish=0.0;
	for(unsigned is=0;is<seqsize;++is) {
		Meff+=(double)(1.0/identity[is]);
                kish+=(double)(1.0/identity[is]) * (double)(1.0/identity[is]);
        }
	double invnorm=1.0/(Meff); 
        kish=Meff*Meff/kish;
	//LEARNING 
	while(t<lenMean+lenObs){
		//METROPOLIS
		int track=0;
		while(track<20){
			unsigned in=0;
			while(in<len){
				iold=fh[in];
				inew=iold;
				while(inew==iold)
					inew=nucl(g);
				Jcum=0.0;	Jcumold=0.0;
				for(unsigned jn=0;jn<len;jn++){
					jold=fh[jn];
					if(in!=jn){
						Jcumold += J[in*len+jn].count[iold][jold];
						Jcum += J[in*len+jn].count[inew][jold];
					}
				}	
				deltaH = h[in].sing[inew]-h[in].sing[iold] + Jcum-Jcumold;
				if(deltaH<0 || exp(-deltaH)>rando(g)){
						fh[in]=inew;
					}
				in++;
			}
		track++;
		}
		t++;
		//CALCOLO PARAMETRI
		if(t<lenMean){	
			for(unsigned in=0;in<len;in++)	
				for(unsigned i=0;i<5;i++){
					if(i==4)h[in].sing[i]=0;
					else
						h[in].sing[i]=h[in].sing[i]+(alfah/std::pow(1+t/(double)tau,exponent))*((fh[in]==i?1:0)-freq[in][i]*invnorm-lh*(h[in].sing[i]));
					for(unsigned jn=0;jn<len;jn++)
						for(unsigned j=0;j<5;j++){
							if(i==4 || j==4)J[in*len+jn].count[i][j]=0;
							else 
							J[in*len+jn].count[i][j]=J[in*len+jn].count[i][j]+(alfaj/std::pow(1+t/(double)tau,exponent))*
							((fh[in]==i?1:0)*(fh[jn]==j?1:0)-histogram[in*len+jn].count[i][j]*invnorm-L*(J[in*len+jn].count[i][j]));
						}		
				}
			//MEDIA PARAMETRI SUI PROCESSORI
			for(unsigned in=0;in<len;in++)
				MPI_Allreduce(&h[in].sing[0],&new_h[in].sing[0],5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
				
			for(unsigned in=0;in<len;in++)for(unsigned jn=0;jn<len;jn++)
				MPI_Allreduce(&J[in*len+jn].count[0][0],&new_J[in*len+jn].count[0][0],25,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);	
			
			for(unsigned in=0;in<len;in++)for(unsigned i=0;i<5;i++){
				h[in].sing[i]=new_h[in].sing[i]/(double)np;
				for(unsigned jn=0;jn<len;jn++)for(unsigned j=0;j<5;j++)
					J[in*len+jn].count[i][j]=new_J[in*len+jn].count[i][j]/(double)np;	
			}
			
			if(t>=lenWait){
				for(unsigned in=0;in<len;in++){
					for(unsigned i=0;i<5;i++){
						h[in].sum[i]+=h[in].sing[i];
						for(unsigned jn=0;jn<len;jn++)
							for(unsigned j=0;j<5;j++)
								J[in*len+jn].mean[i][j]+=J[in*len+jn].count[i][j];
					}
				}
			}	
		}
		//ASSEGNAZIONE VALORI MEDI 
		if((int)t==lenMean){
			for(unsigned in=0;in<len;in++)for(unsigned i=0;i<5;i++){h[in].sing[i]=h[in].sum[i]/(double)(lenMean-lenWait);
				for(unsigned jn=0;jn<len;jn++)for(unsigned j=0;j<5;j++) {J[in*len+jn].count[i][j]=J[in*len+jn].mean[i][j]/(double)(lenMean-lenWait);}}
		}
		//CONTROLLO: RICALCOLO ISOGRAMMI CON PARAMETRI FISSATI E BLOCK ANALYSIS
		if(t>lenMean){
			for(unsigned in=0;in<len;in++)for(unsigned i=0;i<5;i++){
				if(fh[in]==i)
					f[in].sing[i]++;
				for(unsigned jn=0;jn<len;jn++)
					for(unsigned j=0;j<5;j++)
						if(fh[in]==i && fh[jn]==j)
							F[in*len+jn].count[i][j]++;
			}
		}
	}
	//MEDIA SUI PROCESSORI DELLE PROBABILITÀ
	for(unsigned in=0;in<len;in++)for(unsigned jn=0;jn<len;jn++)
		MPI_Allreduce(&F[in*len+jn].count[0][0],&new_F[in*len+jn].count[0][0],25,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	for(unsigned in=0;in<len;in++)
		MPI_Allreduce(&f[in].sing[0],&new_f[in].sing[0],5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	//SCRIVO RISULTATI
	if(!rank){
		std::stringstream fileh;
		fileh<<"convh";
		std::ofstream hh(fileh.str().c_str());
		std::stringstream filej;
		filej<<"convj";
		std::ofstream jj(filej.str().c_str());
		std::stringstream filesing;
		filesing<<"h_val";
		std::ofstream hval(filesing.str().c_str()); 
		std::stringstream filecoup;
		filecoup<<"j_val";
		std::ofstream jval(filecoup.str().c_str()); 
		
		for(unsigned in=0;in<len;in++)
			for(unsigned i=0;i<5;i++){
				f[in].sing[i]=new_f[in].sing[i]/(double)np;
				 hh<<freq[in][i]*invnorm<<" "<<f[in].sing[i]/(double)lenObs-lh*h[in].sing[i]<<"\n"; 
				hval<<h[in].sing[i]<<"\n";	
			}
		jval<<len<<"\n";
		for(unsigned in=0;in<len;in++)	
			for(unsigned jn=0;jn<len;jn++){
				for(unsigned i=0;i<5;i++){
					for(unsigned j=0;j<5;j++){
						F[in*len+jn].count[i][j]=new_F[in*len+jn].count[i][j]/(double)np;
						jj<<histogram[in*len+jn].count[i][j]*invnorm<<" "<<F[in*len+jn].count[i][j]/(double)lenObs-L*J[in*len+jn].count[i][j]<<"\n"; 
						jval<<J[in*len+jn].count[i][j]<<"\n";	
					}
				}
			}
	}
	MPI_Finalize();	
	return 0;
}	
