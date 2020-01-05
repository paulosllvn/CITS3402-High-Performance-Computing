/**
 * CITS3402 Project 1
 * Author : Paul O'Sullivan
 * Student Number: 21492318
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include<sys/time.h>
#include <omp.h>
 
int intFileRead(char* file, int ** nnz, int ** ja,int ** ia, int** row, int * countPtr, int* nonZero){
    int datatypeLine = 0;
    int lengthLine = 1;
    int widthLine = 2;
    char str[512];
    int loop;
    int lengthFloat,widthFloat;
    char* datatype;
 
    fflush(stdout);
 
    FILE *fd= fopen(file, "r");
    if (fd == NULL) {
        printf("Failed to open file\n");
        EXIT_FAILURE;
    }
 
    for(loop = 0;loop<=widthLine;++loop){
        fgets(str,sizeof(str),fd);
 
        if(loop == datatypeLine){
            datatype = str;
 
        }
        if(loop == lengthLine){
            lengthFloat = atoi(str);
 
        }
        if(loop==widthLine){
            widthFloat = atoi(str);
        }
 
    }
 
    int temp;
    *nonZero = 0;
    for(int i = 0; i<lengthFloat;i++){
        for(int j = 0; j<widthFloat; j++){
            fscanf(fd,"%d",&temp);
            if(temp==0){
                continue;
            }
            else{
                (*nonZero)++;
            }
        }
 
    }
    fclose(fd);
    fd = fopen(file, "r");
 
    fgets(str,sizeof(str),fd);
    fgets(str,sizeof(str),fd);
    fgets(str,sizeof(str),fd);
 
    *nnz = malloc(*nonZero * sizeof(int));
    *ia = malloc((lengthFloat+1) * sizeof(int));
    *ja = malloc(*nonZero * sizeof(int));
    *row = malloc(*nonZero * sizeof(int));
   
    int temp2;
    int count = 0;
    int rowCount = 0;
    (*ia)[0]=0;
 
     for(int i = 0; i<lengthFloat;i++){
       
        rowCount = 0;
 
        for(int j = 0; j<widthFloat; j++){
            fscanf(fd,"%d",&temp2);
            if(temp2==0){
            continue;
            }
            else{
                (*nnz)[count] = temp2;
                (*ja)[count] = j;
                (*row)[count]= i; //for COO format
                count++;
                rowCount++;
 
            }
        }
        (*ia)[i+1] = (*ia)[i]+rowCount;
 
    } 
 
    *countPtr = count;
 
    fclose(fd);
 
    return 0;
 
}

int floatFileRead(char* file, float ** nnz_f, int ** ja,int ** ia,int** row,int * countPtr,int * nonZero){
    int datatypeLine = 0;
    int lengthLine = 1;
    int widthLine = 2;
    char str[512];
    int loop;
    int lengthFloat,widthFloat;
    char* datatype;
 
    fflush(stdout);
 
    FILE *fd= fopen(file, "r");
    if (fd == NULL) {
        printf("Failed to open file\n");
        exit(1);
    }
 
    for(loop = 0;loop<=widthLine;++loop){
        fgets(str,sizeof(str),fd);
 
        if(loop == datatypeLine){
            datatype = str;
 
        }
        if(loop == lengthLine){
            lengthFloat = atoi(str);
 
        }
        if(loop==widthLine){
            widthFloat = atoi(str);
        }
 
    }
   
    float temp;
    *nonZero = 0;
    for(int i = 0; i<lengthFloat;i++){
        for(int j = 0; j<widthFloat; j++){
            fscanf(fd,"%f",&temp);
            if(temp==0.0){
                continue;
            }
            else{
                (*nonZero)++;
            }
        }
   
 
        }
   
    fclose(fd);
    fd = fopen(file, "r");
 
    fgets(str,sizeof(str),fd);
    fgets(str,sizeof(str),fd);
    fgets(str,sizeof(str),fd);
 
    *nnz_f = malloc(*nonZero * sizeof(float));
    *ia = malloc((lengthFloat+1) * sizeof(int));
    *ja = malloc(*nonZero * sizeof(int));
    *row = malloc(*nonZero * sizeof(int));
   
    float temp2;
    int count = 0;
    int rowCount = 0;
    (*ia)[0]=0.0;
 
     for(int i = 0; i<lengthFloat;i++){
       
        rowCount = 0;
 
        for(int j = 0; j<widthFloat; j++){
            fscanf(fd,"%f",&temp2);
            if(temp2==0.0){
            continue;
            }
            else{
                (*nnz_f)[count] = temp2;
                (*ja)[count] = j;
                (*row)[count]= i; //for COO format
                count++;
                rowCount++;
 
            }
            }
        (*ia)[i+1] = (*ia)[i]+rowCount;
 
    } 
 
    *countPtr = count;
 
    fclose(fd);
 
    return 0;
 
}

int vectorMultiplyFloat(float* nnzFloatArray, int nonZeroElements, float inputScaler, int numberOfThreads)
{   

    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    #pragma omp parallel for
    for(int i = 0; i<nonZeroElements; i++){
     nnzFloatArray[i] =  nnzFloatArray[i] * inputScaler;
    }
    //To get time of parallelising
    double time = omp_get_wtime()-start;
    
    return 0;
} 

int vectorMultiplyInt(int* nnzIntArray, int nonZeroElements,float inputScaler,int numberOfThreads)
{   
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    #pragma omp parallel for
    for(int i = 0; i<nonZeroElements; i++){
        nnzIntArray[i] = nnzIntArray[i] * inputScaler;
    }
    double time = omp_get_wtime()-start;
    return 0;
}

int additionFloat(float* nnz, float* nnz2,float* nnzAdded,int* ja,int* ja2, int rowLength,int colLength, int nonZero, int nonZero2, int* row, int* row2,int numberOfThreads){

    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    
    #pragma omp parallel for
    for(int i =0;i<rowLength*colLength; i++){
            nnzAdded[i] = 0.0;
    }
    #pragma omp parallel for
    for(int i = 0; i<nonZero; i++){
        int a = row[i]*colLength + ja[i];
        nnzAdded[a] = nnz[i] + nnzAdded[a];
    }
    #pragma omp parallel for
    for(int i = 0; i<nonZero2; i++){
        int a = row2[i]*colLength + ja2[i];
        nnzAdded[a] = nnz2[i] + nnzAdded[a];   
    }

    double time = omp_get_wtime()-start;
    return 2;
    
    }
 
int additionInt(int* nnz, int* nnz2,int* nnzAdded, int* ja,int* ja2, int rowLength,int colLength, int nonZero, int nonZero2, int* row, int* row2,int numberOfThreads){
    
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    
    #pragma omp parallel for
    for(int i =0;i<rowLength*colLength; i++){
            nnzAdded[i] = 0;
    }
    #pragma omp parallel for
    for(int i = 0; i<nonZero; i++){
        int a = row[i]*colLength + ja[i];
        nnzAdded[a] = nnz[i] + nnzAdded[a];
    }
    #pragma omp parallel for
    for(int i = 0; i<nonZero2; i++){
        int a = row2[i]*colLength + ja2[i];
        nnzAdded[a] = nnz2[i] + nnzAdded[a];
        
    }

    double time = omp_get_wtime()-start; 
        return 2;
           
    }
   
 
int traceInt(int* nnz,int* ja,int* ia,int nonZero,int rowLength,int numberOfThreads,int* row){
    
    int sum = 0.0;
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    
    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i<nonZero; i++){
            if(ja[i] == row[i]){
             sum = sum + nnz[i];
            }
        }
    double time = omp_get_wtime()-start;
    // printf("omp time is %g\n",time);
    
    return sum;
}
 
float traceFloat(float* nnz_f,int* ja,int* ia,int nonZero,int rowLength,int numberOfThreads,int* row){
    float sum = 0.0;
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();

    #pragma omp parallel for reduction(+:sum)
    for(int i = 0; i<nonZero; i++){
            if(ja[i] == row[i]){
             sum = sum + nnz_f[i];
            }
        }
    double time = omp_get_wtime()-start;
    // printf("omp time is %g\n",time);
    
    return sum;
}
 
int transposeInt(int* nnz,int* column,int* row,int nonZero,int numberOfThreads){
 
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    #pragma omp parallel for
    for(int i = 0; i<nonZero;i++){
        int temp = column[i];
        column[i] = row[i];
        row[i] = temp;
    }

    double time = omp_get_wtime()-start;
    //printf("omp time is %g\n",time);

    return 1;
  
}
 
int transposeFloat(float* nnz_f,int* column,int* row,int nonZero,int numberOfThreads){
    
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    #pragma omp parallel for
    for(int i = 0; i<nonZero;i++){
        int temp = column[i];
        column[i] = row[i];
        row[i] = temp;
       
    }

    double time = omp_get_wtime()-start;
    //printf("omp time is %g\n",time);

    return 1;
}

int writeoutput(char* filename,char* filename2,int* nnz,float* nnzFloat, int* ia, int* ja,int* temp,float* tempf, int rowLength, int colLength,int nonZero,int type,char* action,int thread,int trace,float tracefloat,double delta,double delta1,int transpose,int* row){
     
    char* studentnumber = "21492318_";
    struct tm *tm;
    time_t t;
    char str_time[100];
    char str_date[100];


    t = time(NULL);
    tm = localtime(&t);

    strftime(str_time, sizeof(str_time), "%H%M_", tm);
    strftime(str_date, sizeof(str_date), "%d%m%Y_", tm);
    
    char combine[1000];
    memset(combine,0,sizeof(combine));
    strcat(combine,studentnumber);
    strcat(combine,str_date);
    strcat(combine,str_time);
    strcat(combine,action);
    strcat(combine,".out");

    FILE* output = fopen(combine,"w");
    fprintf(output,"%s\n",action);
    
    fprintf(output,"%s\n",filename);

    if(filename2 != NULL){
        fprintf(output,"%s\n",filename2);
    }
    fprintf(output,"%d\n",thread);
    if(type==1){
        fprintf(output,"int\n");
    }
    else if(type==2){
        fprintf(output,"float\n");

    }
    fprintf(output,"%d\n",rowLength);
    fprintf(output,"%d\n",colLength);

    
    if(trace != 0){
       fprintf(output,"%d\n",trace); 
    }

    if(tracefloat != 0){
       fprintf(output,"%f\n",tracefloat); 
    }

    if(nnz != NULL){
        if(transpose == 2){
            for(int i = 0; i<colLength*rowLength; i++){
                fprintf(output,"%d ",nnz[i]);
            }
        }
        else if(transpose==1){
            for(int i=0; i<(rowLength); i++){

                for(int j = 0; j<colLength; j++){
                        temp[j] = 0;
                    }
                for(int k = 0; k<nonZero; k++){
                    if(row[k] == i){
                        temp[ja[k]] = nnz[k];
                    }
                }
                for(int g = 0; g<colLength; g++){
                    fprintf(output,"%d ",temp[g]);
                }

            }
        }
        else{
            int k = 0;
            for(int i=0; i<(rowLength); i++){
                int indexVal = (ia)[i+1];
                int indexValPrev = (ia)[i];
                int difference = indexVal-indexValPrev;
                difference = difference + k;
    
    
                for(int j = 0; j<colLength; j++){
                    temp[j] = 0;
                }
        
                for(; k<difference; k++){
                    temp[ja[k]] = nnz[k];
                }
    
                for(int g = 0; g<colLength; g++){
                    fprintf(output,"%d ",temp[g]);

                }
    
            }
        }
    }

    if(nnzFloat != NULL){
        if(transpose == 2){
            for(int i = 0; i<colLength*rowLength; i++){
                if(nnzFloat[i]==0.0){
                        fprintf(output,"%.1f ", nnzFloat[i]);
                    }
                    else{
                        fprintf(output,"%f ",nnzFloat[i]);
                    }
            }
        }
        else if(transpose==1){
            for(int i=0; i<(rowLength); i++){

                for(int j = 0; j<colLength; j++){
                        tempf[j] = 0.0;
                    }
                for(int k = 0; k<nonZero; k++){
                    if(row[k] == i){
                        tempf[ja[k]] = nnzFloat[k];
                    }
                }
                for(int g = 0; g<colLength; g++){
                    if(tempf[g]==0.0){
                        fprintf(output,"%.1f ", tempf[g]);
                    }
                    else{
                        fprintf(output,"%f ",tempf[g]);
                        // printf("printing is %f\n",tempf[g]);
                    }
                }

            }
        }
        else{   
            int k = 0;
            for(int i=0; i<(rowLength); i++){
                int indexVal = (ia)[i+1];
                int indexValPrev = (ia)[i];
                int difference = indexVal-indexValPrev;

                difference = difference+ k;

            
        
                for(int j = 0; j<colLength; j++){
                tempf[j] = 0.0;
                }
            
                while(k<difference){
                    
                    tempf[ja[k]] = nnzFloat[k];
                    k++;
                }


                for(int g = 0; g<colLength; g++){
                    if(tempf[g]==0.0){
                        fprintf(output,"%.1f ", tempf[g]);
                    }
                    else{
                        fprintf(output,"%f ",tempf[g]);
                    }
                }
        
            }
        }
    }  
    fprintf(output,"\n");
    fprintf(output,"%8.6f\n",delta);
    fprintf(output,"%8.6f\n",delta1);
 
        return 0;
}

int writestdout(char* filename,char* filename2,int* nnz,float* nnzFloat, int* ia, int* ja,int* temp,float* tempf, int rowLength, int colLength,int nonZero,int type,char* action,int thread,int trace,float tracefloat,double delta,double delta1,int transpose,int* row){
    char* studentnumber = "21492318_";
    struct tm *tm;
    time_t t;
    char str_time[50];
    char str_date[50];


    t = time(NULL);
    tm = localtime(&t);

    strftime(str_time, sizeof(str_time), "%H%M", tm);
    strftime(str_date, sizeof(str_date), "%d%m%Y_", tm);
    
    char combine[500];
    memset(combine,0,sizeof(combine));
    strcat(combine,studentnumber);
    strcat(combine,str_date);
    strcat(combine,str_time);
    strcat(combine,action);
    strcat(combine,".out");

    printf("%s\n",combine);

    printf("%s\n",action);
    
    printf("%s\n",filename);

    if(filename2 != NULL){

    printf("%s\n",filename2);
    }
    printf("%d\n",thread);
    if(type==1){
        printf("int\n");
    }
    else if(type==2){
        printf("float\n");
    }
    
    printf("%d\n",rowLength);
    printf("%d\n",colLength);
    
    
    if(trace != 0){
       printf("%d\n",trace); 
    }

    if(tracefloat != 0){
       printf("%f\n",tracefloat); 
    }
    
    if(nnz != NULL){
        if(transpose == 2){
            for(int i = 0; i<colLength*rowLength; i++){
                printf("%d ",nnz[i]);
            }
        }
        else if(transpose == 1){
            for(int i=0; i<(rowLength); i++){

                for(int j = 0; j<colLength; j++){
                        temp[j] = 0;
                    }
                for(int k = 0; k<nonZero; k++){
                    if(row[k] == i){
                        temp[ja[k]] = nnz[k];
                    }
                }
                for(int g = 0; g<colLength; g++){
                    printf("%d ",temp[g]);
                }

            }
        }
        else{   
            int k = 0;
            for(int i=0; i<(rowLength); i++){
                int indexVal = (ia)[i+1];
                int indexValPrev = (ia)[i];
                int difference = indexVal-indexValPrev;
                difference = difference + k;
        
        
                for(int j = 0; j<colLength; j++){
                    temp[j] = 0;
                }
            
                for(; k<difference; k++){
                    temp[ja[k]] = nnz[k];
                }
        
                for(int g = 0; g<colLength; g++){
                    printf("%d ",temp[g]);

                }
        
            }
        }
        
    }

    if(nnzFloat != NULL){
        if(transpose == 2){
            for(int i = 0; i<colLength*rowLength; i++){
                if(nnzFloat[i]==0.0){
                        printf("%.1f ", nnzFloat[i]);
                    }
                    else{
                        printf("%f ",nnzFloat[i]);
                    }
            }
        }
        else if(transpose==1){
            for(int i=0; i<(rowLength); i++){

                for(int j = 0; j<colLength; j++){
                        tempf[j] = 0.0;
                    }
                for(int k = 0; k<nonZero; k++){
                    if(row[k] == i){
                        tempf[ja[k]] = nnzFloat[k];
                    }
                }
                for(int g = 0; g<colLength; g++){
                    if(tempf[g]==0.0){
                        printf("%.1f ", tempf[g]);
                    }
                    else{
                        printf("%f ",tempf[g]);
                    }
                }

            }
        }
        else{   
            int k = 0;
            for(int i=0; i<(rowLength); i++){
                int indexVal = (ia)[i+1];
                int indexValPrev = (ia)[i];
                int difference = indexVal-indexValPrev;
                difference = difference + k;
        
        
                for(int j = 0; j<colLength; j++){
                    tempf[j] = 0.0;
                }
            
                for(; k<difference; k++){
                    tempf[ja[k]] = nnzFloat[k];
                }
        
                for(int g = 0; g<colLength; g++){
                    if(tempf[g]==0.0){
                        printf("%.1f ", tempf[g]);
                    }
                    else{
                        printf("%f ",tempf[g]);
                    }
                }
        
            }
        }
    }   
    printf("\n");
    printf("%8.6f\n",delta);
    printf("%8.6f\n",delta1);

        return 0;
}

 
int multiplicationInt(int* nnz,int* nnz2,int* nnzAdded, int* iaAdded, int* jaAdded,int* ja,int* ja2,int* row,int*row2,int rowLength,int rowLength2,int colLength,int colLength2,int nonZero, int nonZero2,int index,int numberOfThreads){
        
        omp_set_num_threads(numberOfThreads);
        double start = omp_get_wtime();
        int sum = 0;
        int indexIA = 0;
        iaAdded[0] = 0;
    
        for(int i = 0; i<rowLength; i++){

            for(int l = 0; l<colLength2; l++){
                
                sum = 0;
                #pragma omp parallel for
                for(int k = 0; k<nonZero;k++){
                    
                    for(int g = 0; g<nonZero2 ;g++ ){
                        if(row[k]==i && ja2[g]==l && ja[k]==row2[g]){
                            sum = sum + nnz[k] * nnz2[g];
                            break;
                        }

                    
                    }
            
                }
                if(sum!=0){
                    indexIA++;
                }

                if(sum != 0){
                    nnzAdded[index] = sum;
                    jaAdded[index] = l;
                    index++;
                }
        
            }
                iaAdded[i+1] = indexIA;  
        }  
            double time = omp_get_wtime()-start;
    
    
    return 0;
 
}

int multiplicationFloat(float* nnzf,float* nnzf2,float* nnzAddedFloat, int* iaAdded, int* jaAdded,int* ja,int* ja2,int* row,int*row2,int rowLength,int rowLength2,int colLength,int colLength2,int nonZero, int nonZero2,int index,int numberOfThreads){

    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();

        float sum = 0.0;
        int indexIA = 0;
        iaAdded[0] = 0;
        
        for(int i = 0; i<rowLength; i++){

            for(int l = 0; l<colLength2; l++){
                
                sum = 0;
                #pragma omp parallel for
                for(int k = 0; k<nonZero;k++){
                    for(int g = 0; g<nonZero2 ;g++ ){
                        if(row[k]==i && ja2[g]==l && ja[k]==row2[g]){
                            
                            sum = sum + nnzf[k] * nnzf2[g];
                            break;
                        }

                    
                    }
            
                }
                if(sum!=0){
                    indexIA++;
                }

                if(sum != 0){
                    nnzAddedFloat[index] = sum;
                    jaAdded[index] = l;
                    index++;
                }
        
            }
                iaAdded[i+1] = indexIA;  
        }  
        double time = omp_get_wtime()-start;
    return 0;
 
} 


int getType(char* file){
    char str[512];
    FILE *fp= fopen(file, "r");
    fgets(str,sizeof(str),fp);
    fclose(fp);
    if(!strncmp(str,"int",3)){
        return 1;
    }
    else{
        return 2;
    }
 
}
 
int getColRow(char* file, int* rowLength, int* colLength){
    char str[512];
    FILE *fp= fopen(file, "r");
    for(int i = 0; i<3; i++){
        fgets(str,sizeof(str),fp);
        if(i == 1){
            (*rowLength) = atoi(str);
    
        }
        if(i==2){
            (*colLength) = atoi(str);
        }  
    }
 return 0;
}
 
int main(int argc, char **argv) {
 
    char *action = NULL;
    char *filename1=NULL;
    char* filename2 = NULL;
    bool Log = false;
    float scaler = 1.0;
    int threads = 1;
    int * nnz = NULL;
    float * nnz_f = NULL;
    int * ja = NULL;
    int * ia = NULL;
    int* row = NULL;
    int* row2 = NULL;
    int * nnz2 = NULL;
    float * nnz_f2 = NULL;
    int * ja2 = NULL;
    int * ia2 = NULL;
    int loopCount,loopcount2,nonZero,nonZero2,colLength,colLength2,rowLength,rowLength2,trace = 0;
    int transpose = 0;
    int mainIndex =0;
    int nnzAddedSize = 0;
    float tracefloat = 0.0;
    int* iaTranspose = NULL;
    int* jaTranspose = NULL;
    int* nnzTranspose = NULL;
    int* nnzAdded = NULL;
    int* iaAdded = NULL;
    int* jaAdded = NULL;
    float* nnzAddedFloat = NULL;
    int* temp = NULL;
    float* tempf = NULL;
   
    int bigNum = 100000;
    double delta = 0;
    double delta1 = 0;
    int b[bigNum];
    int i,sum =0;


    struct timeval start, end;
    gettimeofday(&start, NULL);

    for( i = 0; i < bigNum; i++)
    {
        sum += b[i];
    }
 
 
    // Reads the comination of arguments skipping the first
    for (int i = 1; i < argc; i++){
 
        if((strcmp(argv[i], "--sm") ==0 || strcmp(argv[i], "--tr" ) ==0 || strcmp(argv[i], "--ad" ) ==0 || strcmp(argv[i], "--ts" ) == 0|| strcmp(argv[i], "--mm" )==0)){
            action = argv[i]+2;
            }
    
        if(strcmp(argv[i],"-l")==0){
            Log = true;
        }
    
        if(strcmp(argv[i],"-f")==0){
            if(argv[i+1]){
                filename1 = argv[i+1];
            }
        }
    
        if((strcmp(action,"ad")==0 || strcmp(action,"mm")==0) && strcmp(argv[i],"-f")==0){
            filename2 = argv[i+2];
        }
        
    
        if(strcmp(argv[i],"--sm")==0){
            if(argv[i+1]){
                scaler = atof(argv[i+1]);
                }
        }
    
        if(strcmp(argv[i],"-t")==0) {
            if(argv[i+1]){
                threads = atoi(argv[i+1]);
            }
        }
 
    }

    int a = getType(filename1);

    
 
    if (strcmp(action, "ad")==0){
        if(a==1){
            intFileRead(filename1, &nnz, &ja, &ia,&row, &loopCount,&nonZero);
            intFileRead(filename2, &nnz2, &ja2, &ia2,&row2, &loopcount2,&nonZero2);
            
            //Get time to read in matrix
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;


            getColRow(filename1,&rowLength,&colLength);
            getColRow(filename1,&rowLength2,&colLength2);
            
            
            if(rowLength!=rowLength2 || colLength!=colLength2){
                printf("ERROR, matrices must be the same dimensions");
                    EXIT_FAILURE;
            } 
            else{
                nnzAdded = malloc((rowLength*colLength) * (sizeof(int)));
                iaAdded = malloc((rowLength+1) * (sizeof(int)));
                jaAdded = malloc((nonZero+nonZero2) * (sizeof(int)));
                temp = malloc((colLength) * (sizeof(int)));
                 gettimeofday(&start,NULL);
                int addition = additionInt(nnz,nnz2,nnzAdded,ja,ja2,rowLength,colLength,nonZero,nonZero2,row,row2,threads);
                
                //Get time to perform operation on matrix
                gettimeofday(&end, NULL);
                delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u +
                    end.tv_usec - start.tv_usec) / 1.e6;

                if(Log == true){
    
                    writeoutput(filename1,filename2,nnzAdded,nnz_f,iaAdded,jaAdded,temp,tempf,rowLength,colLength,mainIndex,a,action,threads,trace,tracefloat,delta,delta1,addition,row); 
                
                }
                else{
                
                    writestdout(filename1,filename2,nnzAdded,nnz_f,iaAdded,jaAdded,temp,tempf,rowLength,colLength,mainIndex,a,action,threads,trace,tracefloat,delta,delta1,addition,row);
                }
                free(temp);      
            }  
                
           
        }
        
        if(a==2){
            floatFileRead(filename1, &nnz_f, &ja, &ia,&row, &loopCount,&nonZero);
            floatFileRead(filename2, &nnz_f2, &ja2, &ia2,&row2, &loopcount2,&nonZero2);
            
            //Get time to read in matrix
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

            getColRow(filename1,&rowLength,&colLength);
            getColRow(filename2,&rowLength2,&colLength2);
            
            if(rowLength!=rowLength2 || colLength!=colLength2){
                printf("ERROR, matrices must be the same dimensions");
                    EXIT_FAILURE;
            } 
            else{
                nnzAddedFloat = malloc((rowLength*colLength) * (sizeof(float)));
                iaAdded = malloc((rowLength+1) * (sizeof(int)));
                jaAdded = malloc((nonZero+nonZero2) * (sizeof(int)));
                tempf = malloc((colLength) * (sizeof(float)));
                gettimeofday(&start,NULL);
                int add = additionFloat(nnz_f,nnz_f2,nnzAddedFloat,ja,ja2,rowLength,colLength,nonZero,nonZero2,row,row2,threads);
                
                gettimeofday(&end, NULL);
                delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

                
                if(Log == true){ 
                    writeoutput(filename1,filename2,nnz,nnzAddedFloat,iaAdded,jaAdded,temp,tempf,rowLength,colLength,mainIndex,a,action,threads,trace,tracefloat,delta,delta1,add,row);
                }
                else{ 
                    writestdout(filename1,filename2,nnz,nnzAddedFloat,iaAdded,jaAdded,temp,tempf,rowLength,colLength,mainIndex,a,action,threads,trace,tracefloat,delta,delta1,add,row);
                }
                free(tempf);
            }
        }
    
        free(nnz_f); free(nnz_f2); free(ia);free(ia2); free(ja); free(ja2); free(nnzAddedFloat); free(iaAdded); free(jaAdded);free(row);
        free(row2);  free(nnzAdded); free(nnz); free(nnz2); 
                
            
    }
    
    
   
    if (strcmp(action, "ts")==0)
    {
        if (a==1){
            intFileRead(filename1, &nnz, &ja, &ia,&row, &loopCount,&nonZero);

            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

            getColRow(filename1,&rowLength,&colLength);
            gettimeofday(&start,NULL);

            int trans = transposeInt(nnz,ja,row,nonZero,threads);
            
            gettimeofday(&end, NULL);
            delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

            temp = malloc((rowLength) * (sizeof(int)));
            if(Log==true){
                writeoutput(filename1,filename2,nnz,nnz_f,iaTranspose,ja,temp,tempf,colLength,rowLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,trans,row);
            }
            else{
                writestdout(filename1,filename2,nnz,nnz_f,iaTranspose,ja,temp,tempf,colLength,rowLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,trans,row);
            }
            free(nnz); free(ia); free(ja); free(row);
            free(temp);
 
       
        }
        if(a==2){
            floatFileRead(filename1, &nnz_f, &ja, &ia,&row, &loopCount,&nonZero);
            
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            getColRow(filename1,&rowLength,&colLength);
            gettimeofday(&start,NULL);

            int trans = transposeFloat(nnz_f,ja,row,nonZero,threads);
            
            gettimeofday(&end, NULL);
            delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

            tempf = malloc((colLength) * (sizeof(float)));
            if(Log==true){
                writeoutput(filename1,filename2,nnzTranspose,nnz_f,iaTranspose,ja,temp,tempf,colLength,rowLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,trans,row);
            }
            else{
                writestdout(filename1,filename2,nnzAdded,nnz_f,iaTranspose,ja,temp,tempf,colLength,rowLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,trans,row);
            }
            free(nnz_f); free(ia); free(ja); free(row);
            free(tempf);
        }
             
    }
 
   
    if(strcmp(action, "tr")==0){
        if(a==1){
            intFileRead(filename1, &nnz, &ja, &ia,&row, &loopCount,&nonZero);
            
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            getColRow(filename1,&rowLength,&colLength);
            gettimeofday(&start,NULL);
            trace = traceInt(nnz,ja,ia,nonZero,rowLength,threads,row);
            
            
            gettimeofday(&end, NULL);
            delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            temp = malloc((colLength) * (sizeof(int)));
            if(Log==true){
                writeoutput(filename1,filename2,nnzTranspose,nnz_f,iaTranspose,jaTranspose,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
            else{
                writestdout(filename1,filename2,nnzAdded,nnz_f,iaAdded,jaAdded,temp,tempf,rowLength,colLength,mainIndex,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
            free(nnz); free(ia); free(ja); free(row);
            free(temp);
        }

        else if(a==2){
            floatFileRead(filename1, &nnz_f, &ja, &ia,&row, &loopCount,&nonZero);
            
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            getColRow(filename1,&rowLength,&colLength);
            gettimeofday(&start,NULL);
            tracefloat = traceFloat(nnz_f,ja,ia,nonZero,rowLength,threads,row);
            
            gettimeofday(&end, NULL);
            delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            tempf = malloc((colLength) * (sizeof(float)));
            if(Log==true){
                writeoutput(filename1,filename2,nnzTranspose,nnz_f2,iaTranspose,jaTranspose,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
            else{
                writestdout(filename1,filename2,nnzTranspose,nnz_f2,iaTranspose,jaTranspose,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
             free(nnz_f); free(ia); free(ja); free(row);
             free(tempf);

        }
    }
 
    if(strcmp(action,"sm")==0){
        if(a == 1){
            intFileRead(filename1, &nnz, &ja, &ia,&row, &loopCount,&nonZero);
            
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            getColRow(filename1,&rowLength,&colLength);
            gettimeofday(&start,NULL);
            vectorMultiplyInt(nnz,nonZero,scaler,threads);
            
            gettimeofday(&end, NULL);
            delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            temp = malloc((colLength) * (sizeof(int)));
            if(Log==true){ 
                writeoutput(filename1,filename2,nnz,nnz_f,ia,ja,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
            else{
                writestdout(filename1,filename2,nnz,nnz_f,ia,ja,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
            free(nnz); free(ia); free(ja); free(row); free(temp); free(row);

           
        }
        else if(a==2){
 
            floatFileRead(filename1, &nnz_f, &ja, &ia,&row, &loopCount,&nonZero);
            
            gettimeofday(&end, NULL);
            delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            getColRow(filename1,&rowLength,&colLength);
            gettimeofday(&start,NULL);
            vectorMultiplyFloat(nnz_f,nonZero,scaler,threads);
            
            gettimeofday(&end, NULL);
            delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            
            tempf = malloc((colLength) * (sizeof(float)));
            if(Log==true){
                writeoutput(filename1,filename2,nnz,nnz_f,ia,ja,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }  
            else{
                writestdout(filename1,filename2,nnz,nnz_f,ia,ja,temp,tempf,rowLength,colLength,nonZero,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
            }
            free(nnz_f); free(ia); free(ja); free(row); free(tempf); free(row);
        }
           
    }
    if(strcmp(action,"mm")==0){
            if(a==1){
                intFileRead(filename1, &nnz, &ja, &ia,&row, &loopCount,&nonZero);
                intFileRead(filename2, &nnz2, &ja2, &ia2,&row2, &loopCount,&nonZero2);
                
                gettimeofday(&end, NULL);
                delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                
                getColRow(filename1,&rowLength,&colLength);
                getColRow(filename2,&rowLength2,&colLength2);

                if(colLength != rowLength2){
                    printf("ERROR, unequal row and column length");
                    EXIT_FAILURE;
                }
                else{
                    nnzAdded = malloc((colLength*rowLength2) * sizeof(int));
                    jaAdded = malloc((colLength*rowLength2) * sizeof(int));
                    iaAdded = malloc((colLength+1) * sizeof(int));
                    gettimeofday(&start,NULL);
                    multiplicationInt(nnz,nnz2,nnzAdded,iaAdded,jaAdded,ja,ja2,row,row2,rowLength,rowLength2,colLength,colLength2,nonZero,nonZero2,nnzAddedSize,threads);
                    
                    gettimeofday(&end, NULL);
                    delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

                    temp = malloc((colLength+1) * (sizeof(int)));
                    if(Log==true){
                        writeoutput(filename1,filename2,nnzAdded,nnzAddedFloat,iaAdded,jaAdded,temp,tempf,rowLength2,colLength,nnzAddedSize,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
                    }
                    else{
                        writestdout(filename1,filename2,nnzAdded,nnzAddedFloat,iaAdded,jaAdded,temp,tempf,rowLength2,colLength,nnzAddedSize,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
                    }
                    
                    free(nnz); free(ia); free(ja); free(row); free(temp); free(row2); free(nnz2); free(ja2); free(ia2);
                    free(nnzAdded); free(jaAdded); free(iaAdded);

                   
                }
            }
            if(a==2){
                floatFileRead(filename1, &nnz_f, &ja,&ia,&row,&loopCount,&nonZero);
                floatFileRead(filename2, &nnz_f2, &ja2, &ia2,&row2, &loopCount,&nonZero2);
                
                gettimeofday(&end, NULL);
                delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                
                getColRow(filename1,&rowLength,&colLength);
                getColRow(filename2,&rowLength2,&colLength2);
                
                if(colLength != rowLength2){
                    printf("ERROR, unequal row and column length");
                    EXIT_FAILURE;
                }
                else{
                    nnzAddedFloat = malloc((colLength*rowLength2) * sizeof(float));
                    jaAdded = malloc((colLength*rowLength2) * sizeof(int));
                    iaAdded = malloc((colLength+1) * sizeof(int));
                    gettimeofday(&start,NULL);
                    multiplicationFloat(nnz_f,nnz_f2,nnzAddedFloat,iaAdded,jaAdded,ja,ja2,row,row2,rowLength,rowLength2,colLength,colLength2,nonZero,nonZero2,nnzAddedSize,threads);
                    
                    gettimeofday(&end, NULL);
                    delta1 = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
                
                    tempf = malloc((colLength+1) * (sizeof(float)));
                    
                    if(Log==true){
                        writeoutput(filename1,filename2,nnzAdded,nnzAddedFloat,iaAdded,jaAdded,temp,tempf,rowLength2,colLength,nnzAddedSize,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
                    }
                    else{
                        writestdout(filename1,filename2,nnzAdded,nnzAddedFloat,iaAdded,jaAdded,temp,tempf,rowLength2,colLength,nnzAddedSize,a,action,threads,trace,tracefloat,delta,delta1,transpose,row);
                    } 
                    free(nnz_f); free(ia); free(ja); free(row); free(tempf); free(row2); free(nnz_f2); free(ja2); free(ia2);
                    free(nnzAddedFloat); free(jaAdded); free(iaAdded);  
                }
                
            } 
    }
}
