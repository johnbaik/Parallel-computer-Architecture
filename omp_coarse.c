#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define TRUE 1
#define FALSE 0



int** initialize(char seq1[],char seq2 []);
int** ScoreTable(char seq1[],char seq2 [],int**M);
int Traceback(char seq1[],char seq2 [],int** M);


char*name;
char*path;
int threads;

int GAP;
int MATCH;
int MISSMATCH;
int SCORE;
FILE *fw;



long count_cells=0;//
long traceback_steps=0;//
double total_exec_time=0;//
double cells_calc_time=0;
double traceback_time=0;//
double CPUS_TOTAL=0;//
double CPUS_cell_calc=0;//





int cont = TRUE;

char dash = '-';
int check;


double gettime(void)
{
struct timeval ttime;
gettimeofday(&ttime, NULL);
return ttime.tv_sec+ttime.tv_usec * 0.000001;
}


int Align(int PosA,int PosB,char seq1[],char seq2 [],int**M,char tracebackQ[],char tracebackD[] ) {
    int ii,jj,MaxAcounter=1,MaxBcounter=1;
	/*Function Variables*/
	int relmax = -1;		//hold highest value in sub columns and rows
	int relmaxpos[2];		//holds position of relmax
        double s =gettime();
	relmaxpos[0]=PosA;
        relmaxpos[1]=PosB;
        double st =gettime();
        cells_calc_time=cells_calc_time+(st-s);
        count_cells++;count_cells++;
        
        
	while(1) {	//until the diagonal of the current cell has a value of zero
		PosA=relmaxpos[0];
                PosB=relmaxpos[1];
                relmax=-1;
                
                //traceback_steps++;
                if(M[PosA-1][PosB-1] == 0) {
                        cont = FALSE;
                        //startD=PosB;
                        break;
                    }
		/*Find relmax in sub columns and rows*/
		for(ii=PosA; ii>0; --ii) {
	
			if(relmax < M[ii-1][PosB-1]) {
                                double s =gettime();
                                relmax = M[ii-1][PosB-1];
                                relmaxpos[0]=ii-1;
                                relmaxpos[1]=PosB-1;
                                count_cells++;count_cells++;
                                double st =gettime();
                                cells_calc_time=cells_calc_time+(st-s);
			}
		}

		for(jj=PosB; jj>0; --jj) {

			if(relmax < M[PosA-1][jj-1]) {
                                double s =gettime();
                                relmax = M[PosA-1][jj-1];
                                relmaxpos[0]=PosA-1;
                                relmaxpos[1]=jj-1;
                                double st =gettime();
                                cells_calc_time=cells_calc_time+(st-s);
                                count_cells++;count_cells++;
			}
		}

		/*Align strings to relmax*/
		if((relmaxpos[0] == PosA-1) && (relmaxpos[1] == PosB-1)) {	//if relmax position is diagonal from current position simply align and increment counters
                        double s =gettime();
			tracebackQ[MaxAcounter] = seq1[relmaxpos[0]-1];
			++MaxAcounter;
			tracebackD[MaxBcounter] = seq2[relmaxpos[1]-1];
			++MaxBcounter;
                         double st =gettime();
                         cells_calc_time=cells_calc_time+(st-s);
                          count_cells++;count_cells++;
		}

		else {

			if((relmaxpos[1] == PosB-1) && (relmaxpos[0] != PosA-1))
                        {	//maxB needs at least one '-'
				
                                for(ii=PosA-1; ii>relmaxpos[0]-1; --ii) {	//for all elements of strA between PosA and relmaxpos[0]
                                        double s =gettime();
                                        tracebackQ[MaxAcounter]= seq1[ii-1];
                                        double st =gettime();
                                        cells_calc_time=cells_calc_time+(st-s);
                                                                count_cells++;
                                        ++MaxAcounter;
                                }

                                for(jj=PosA-1; jj>relmaxpos[0]; --jj) {	//set dashes to MaxB up to relmax
                                         double s =gettime();
                                        tracebackD[MaxBcounter] = dash;
                                         double st =gettime();
                                        cells_calc_time=cells_calc_time+(st-s);
                                        count_cells++;
                                        ++MaxBcounter;
                                }

                                double s =gettime();
                                tracebackD[MaxBcounter] = seq2[relmaxpos[1]-1];	//at relmax set pertinent strB value to MaxB
                                double st =gettime();
                                cells_calc_time=cells_calc_time+(st-s);
                                count_cells++;
                                ++MaxBcounter;
                                
			}

			if((relmaxpos[0] == PosA-1) && (relmaxpos[1] != PosB-1)) {	//MaxA needs at least one '-'

                                    for(jj=PosB-1; jj>relmaxpos[1]-1; --jj) {	//for all elements of strB between PosB and relmaxpos[1]
                                            double s =gettime();
                                            tracebackD[MaxBcounter] = seq2[jj-1];
                                            double st =gettime();
                                            cells_calc_time=cells_calc_time+(st-s);
                                            count_cells++;
                                            ++MaxBcounter;
                                    }
                                    for(ii=PosB-1; ii>relmaxpos[1]; --ii) {		//set dashes to strA
                                             double s =gettime();
                                            tracebackQ[MaxAcounter] = dash;
                                            double st =gettime();
                                            cells_calc_time=cells_calc_time+(st-s);
                                            count_cells++;
                                            ++MaxAcounter;
                                    }
                                    double s =gettime();
                                    tracebackQ[MaxAcounter] = seq1[relmaxpos[0]-1];
                                    double st =gettime();
                                    cells_calc_time=cells_calc_time+(st-s);
                                    count_cells++;
                                    ++MaxAcounter;
			}
		}

		
	}
               traceback_steps=traceback_steps+MaxAcounter;
	return(relmaxpos[1]);
}


int main(int argc, char**args){
   
    name=args[1];
    path=args[2];
    MATCH=atoi(args[3]);    
    MISSMATCH=atoi(args[4]);
    GAP=atoi(args[5]);
    threads=atoi(args[6]);
    
    int flag=0;
    int cc=0;
    char tmp[100]; //First Sequence
    char tmp1[100]; //First Sequence
    
char *seq1;//[100000]; //First Sequence
char *seq2;//[100000]; //First Sequence
    
    
    double start_time =gettime();

    
seq1 = (char *) malloc(1);
       seq2 = (char *) malloc(1);
       
       seq1[0]='\0';
       seq2[0]='\0';


    
    
FILE *fp=fopen( path, "r" );

fw = fopen(name, "w");
if (fp==NULL || fw==NULL){
    return -1;}



fscanf(fp,"%s",tmp);

int pairs=0;
fscanf(fp,"%d",&pairs);


char **Q;
Q = (char **) malloc(pairs*sizeof(char*));
char **D;
D = (char **) malloc(pairs*sizeof(char*));


    while(1){
        fscanf(fp,"%s",tmp);
        if (tmp[0]=='Q' && tmp[1]==':'){

            break;
        }


    }



int k=0;

while (flag!=-1){
   
                        while(1){
                            if(fscanf(fp,"%s",tmp1)!=1){
                                flag=-1;
                            break;}

                            if (tmp1[0]=='D' && tmp1[1]==':'){
                                break;
                            }
                            else {
                                double s =gettime();
                                seq1 =  (char*)realloc(seq1, strlen(tmp1)+strlen(seq1)+1);

                                strcat(seq1,tmp1);
                                double st =gettime();
                                cells_calc_time=cells_calc_time+(st-s);

                            }
                        }
            
                        Q[k]=seq1;
         

                    while(1){
                            if(fscanf(fp,"%s",tmp1)!=1){
                                flag=-1;
                                break;}
                            if (tmp1[0]=='Q' && tmp1[1]==':'){

                                break;                
                            }
                            else { 
                                    double s =gettime();

                                                seq2 =  (char*)realloc(seq2, strlen(tmp1)+strlen(seq2)+1);
                                                
                                                 strcat(seq2,tmp1);
                                    double st =gettime();
                                    cells_calc_time=cells_calc_time+(st-s);

                                }
                         }

           D[k]=seq2;

                 seq1=NULL;
                seq2=NULL;
                seq1 = (char *) malloc(1);
               seq2 = (char *) malloc(1);

               seq1[0]='\0';
               seq2[0]='\0';       
            
 k=k+1;

     
  }// -----------------telos while(1) -----------------------

int i;
int** m;
 int** s;

 
 #pragma omp parallel for private(i,m,s) num_threads(threads)

            for( i = 0; i <k ; i++){
               
                m=initialize(Q[i],D[i]);            
                 s= ScoreTable(Q[i],D[i],m);
                #pragma omp critical
                {
                 Traceback(Q[i],D[i],s);
                }
                
                 
            }


               count_cells=count_cells+strlen(seq1)+strlen(seq2);
           

                free(seq1);
                free(seq2);

                seq1=NULL;
                seq2=NULL;
                seq1 = (char *) malloc(1);
               seq2 = (char *) malloc(1);

               seq1[0]='\0';
               seq2[0]='\0';



    double end_time =gettime();
    total_exec_time=end_time-start_time;
    pairs=k;
      printf("\n-------- Results  ----------");
  
    printf("\nA: pairs= %d",pairs);
    printf("\nB: total_cells = %ld",count_cells);    
    printf("\nC: total traceback steps = %ld",traceback_steps);
    printf("\nD: total exec_time = %f",total_exec_time);
    printf("\nE: total cells calc time = %f",cells_calc_time);
    printf("\nF: total traceback time = %f",traceback_time);
    CPUS_TOTAL=count_cells/total_exec_time;// KELIA POU PERNOUN TIMH STHN MONADA TOU XRONOU
    printf("\nG: CPUS TOTAL TIME = %f",CPUS_TOTAL);
    CPUS_cell_calc=count_cells/cells_calc_time;
    printf("\nH: CPUS cell calc = %f ",CPUS_cell_calc);

    printf("\n----------------------------\nProgram teminated...");



    fclose(fp);   
    //free(fp);
    
    fclose(fw);   
    //free(fw);
    
    return 0;
}

//Step 1: Initialize the table
int** initialize(char seq1[],char seq2 []){
    int **M;
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);
    int val = 0;

    M = malloc((seq1len+1) * sizeof(int*)); // +1 giati vazei midenika stin 1h gramh kai stilh
	if(M == NULL)
		{
		fprintf(stderr, "out of memory\n");
		return NULL;
		}
	for(int i = 0; i <= seq1len; i++)
		{
		M[i] = malloc((seq2len+1) * sizeof(int));
		if(M[i] == NULL)
			{
			fprintf(stderr, "out of memory\n");
			return NULL;
			}
		}
    

	for(int i = 0; i <= seq1len; i++)
		{
		for(int j = 0; j <= seq2len; j++)
			M[i][j] = 0;
		}
	

    return M;
}

                    
//Step 2: ScoreTable creation
int** ScoreTable(char seq1[],char seq2 [],int **M){
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);
                        double s =gettime();
    for (int i = 1; i <= seq1len ; i++)
    {
        for (int j = 1; j <= seq2len ; j++)
        {
            int scoreDiag = 0;
            int scoreLeft=0;
            int scoreUp=0;
            int maxScore=0;
            
            count_cells++;//-------------------------------------------------------------------------------------------------------
            
            if (seq1[i - 1] == seq2[j - 1]){
                scoreDiag = M[i - 1][j - 1] + MATCH;
            }
            else{
                scoreDiag = M[i - 1][j - 1] + MISSMATCH;
            }

            scoreLeft = M[i][j - 1] + GAP;
            scoreUp =  M[i - 1][j] + GAP;

            maxScore = MAX(MAX(scoreDiag, scoreLeft), scoreUp);

            if(maxScore <= 0){
                 M[i][j] = 0;
            }
            else{
                M[i][j] = maxScore;
                

            }
        }
    }
                    double st =gettime();
                    
                    cells_calc_time=cells_calc_time+(st-s);
                
                    return M;
}

//Step 3: PrintTable function
int Traceback(char seq1[],char seq2 [],int** M){
   
    int maxScore=0;
    int **maxScores;
char tracebackD[100700];
char tracebackQ[100700];
    char tmp;
    int count=0;
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);
    
    maxScores = malloc(2 * sizeof(int *));
    
 int ii=0,jj=0;   
    
    
    int MaxAcounter = 1;	//MaxA counter	
    int MaxBcounter = 1;	//MaxB counter
    
        
    
    
    for (int i = 0; i <= seq1len ; i++){
        for (int j = 0; j <= seq2len ; j++){

            if (maxScore<=M[i][j]) {
                
                if(maxScore==M[i][j]){
                                    double s =gettime();
                        maxScores = (int **) realloc(maxScores, (count +1) * sizeof(*maxScores));

                            if(maxScores == NULL)
                            {
                                fprintf(stderr, "out of memory\n");
                                return -1;
                            }

                        maxScores[count] = (int *)malloc(2 * sizeof(int));
                        (maxScores[count])[0]=i;
                        (maxScores[count])[1]=j;
                    double st =gettime();
                    cells_calc_time=cells_calc_time+(st-s);
                          count_cells++;
                           count++;       
                }
                else{
                    
                      for(int i = 0; i < count; i++){
                            free(maxScores[i]);
                            count_cells--;

                      }
                                        double s =gettime();
                      count=0;
                      maxScores = (int **) realloc(maxScores, (count + 1) * sizeof(*maxScores));
                      maxScores[count] = (int *)malloc(2 * sizeof(int));
                      (maxScores[count])[0]=i;
                      (maxScores[count])[1]=j;
                        
                      maxScore=M[i][j];
                    double st =gettime();
                    cells_calc_time=cells_calc_time+(st-s);


                     count=1; 
                      
                
                }
            
            }

        }

    }

     
     
//----------------------- trace back-----------------------

  double tmptime1=gettime();   
 
     char OptA[100700];
     char OptB[100700];
     int startD=0;
          tracebackQ[0] = '\0';
	tracebackD[0]  = '\0';
         
     fprintf(fw,"%s %s\n%s %s\n","Q:",seq1,"D:",seq2);

     for(int k=0; k<count; k++){
 
        OptA[0] = '\0';
	OptB[0]  = '\0';
         int K=0;
                 double s =gettime();
	tracebackQ[0] = seq1[(maxScores[k])[0]-1];
	tracebackD[0] = seq2[(maxScores[k])[1]-1];
                 double st =gettime();
                    cells_calc_time=cells_calc_time+(st-s); 
        count_cells++;
        count_cells++;
        
         startD= Align((maxScores[k])[0],(maxScores[k])[1],seq1,seq2,M,tracebackQ,tracebackD);
  ;
                             double ss =gettime();
         for(ii = strlen(tracebackQ)-1; ii >= 0; --ii) {
		OptA[K] = tracebackQ[ii];
                count_cells++;//------------------------------------------------------------------------
		++K;
	}
                    double sts =gettime();
                    cells_calc_time=cells_calc_time+(st-s);

	K=0;
                            double sss =gettime();
	for(jj=strlen(tracebackD)-1; jj >= 0; --jj) {
		OptB[K] = tracebackD[jj];
                count_cells++;
		++K;
	}
                    double stss =gettime();
                    cells_calc_time=cells_calc_time+(st-s);
         //printf("optB: %s\n\n",OptB);
        
        fprintf(fw,"MATCH %d [ Score: %d, Start:%d , Stop: %d]\nQ: %s\nD: %s\n\n",k+1,maxScore,startD,(maxScores[k])[1],OptA,OptB);

     }    
 
double tmptime2=gettime();  
traceback_time=traceback_time+(tmptime2-tmptime1);
     
     
   //------------------------------------------------------------------
 
    
     for(int i = 0; i <= strlen(seq1); i++)
                        free(M[i]);

for(int i = 0; i < count; i++){
            free(maxScores[i]);}
                   
                free(maxScores);
                free(M);
     return count;
}
