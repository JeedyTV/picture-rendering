#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

enum coord {x,y,z};
enum point {a,b,c};

typedef struct pixel_ {
    unsigned char R;
    unsigned char G;
    unsigned char B;
    float depth;
}Pixel;

// Function that return
// dot product of two vector array.
float dotProduct(float vect_A[], float vect_B[]);

// Function to find
// cross product of two vector array.
void crossProduct(float vect_A[], float vect_B[], float cross_P[]);

//Function that retrun the matrix product.
void mulMat(int m_1,int n_1,float mat1[m_1][n_1],int m_2,int n_2,float mat2[m_2][n_2],float rslt[m_1][n_2]);

//Return the norm of a vector
float norm(float vect[]);

//return the vector normal of a triangle
void getNormal(float A[],float B[],float C[],float normal[]);

int main(int argc, char * argv []){
    double start = omp_get_wtime(); //counter for total time of computation
    double endLoop; //counter to parralel time of computation
    double startLoop; //counter to parralel time of computation
    
    assert(argc == 3);
    
    unsigned char CoR;
    unsigned char CoG;
    unsigned char CoB;

    unsigned char CbgR;
    unsigned char CbgG;
    unsigned char CbgB;

    float pMx;
    float pMy;
    float pMz;

    unsigned int h;
    unsigned int w;

    float theta_y;
    float n;
    float f;


    float eX[] = {1,0,0};
    float eY[] = {0,1,0};
    float eZ[] = {0,0,1};
    
    FILE * param = fopen(argv[2], "r");
    assert(param != NULL);
    fscanf(param, "%hhu %hhu %hhu\n%hhu %hhu %hhu\n%f %f %f\n%u %u\n%f %f %f",&CoR,&CoG,&CoB,&CbgR,&CbgG,&CbgB,&pMx,&pMy,&pMz,&h,&w,&theta_y,&n,&f );
    fclose(param);


    float D = sqrt(pMx*pMx+pMy*pMy+pMz*pMz);    
    float Dg = sqrt(pMx*pMx+pMz*pMz);

    float eU[3];
    for(int i=0;i<3;i++) eU[i] = (1/Dg)*(pMz*eX[i]-pMx*eZ[i]);
    float eV[3]; 
    for(int i=0;i<3;i++) eV[i] = (1/(Dg*D))*(-pMx*pMy*eX[i]+(pMz*pMz+pMx*pMx)*eY[i]-pMy*pMz*eZ[i]);
    float eW[3];
    for(int i=0;i<3;i++) eW[i] = (1/D)*(pMx*eX[i]+pMy*eY[i]+pMz*eZ[i]);
    
    float L[] = {-1/sqrt(3),-1/sqrt(3),-1/sqrt(3)}; 

    float OM[] = {pMx,pMy,pMz};

    float r = (float)w/(float)h;
    float Sy = 1/tan(theta_y/2);
    float Sx = Sy/r;


    FILE *ptr = fopen(argv[1],"rb");
    assert(ptr != NULL);
    unsigned T;
    fread(&T,sizeof(T),1,ptr);
    float ***model= malloc(T*sizeof(float**));
    assert(model != NULL);
    for(int i =0;i<T;i++){
        model[i] = malloc(3*sizeof(float*));
        assert(model[i] != NULL);
    }
    for(int i =0;i<T;i++){
        for(int j=0;j<3;j++){
            model[i][j] = malloc(3*sizeof(float));
            assert(model[i][j] != NULL);
        } 
    }

    float buffer[9];

    for(int count = 0; fread(buffer,sizeof(float),9,ptr) > 0;count++){
        
        model[count][0][0]=buffer[0];
        model[count][0][1]=buffer[1];
        model[count][0][2]=buffer[2];
        model[count][1][0]=buffer[3];
        model[count][1][1]=buffer[4];
        model[count][1][2]=buffer[5];
        model[count][2][0]=buffer[6];
        model[count][2][1]=buffer[7];
        model[count][2][2]=buffer[8];
    }
    
    fclose(ptr);
    
//allocate an array of T for the shade
    float *shade = malloc(sizeof(float)*T);
    assert(shade != NULL);
    float normal[3]; // modif stef
    for(int i =0;i<T;i++){
        
        // compute normal of the triangle
        getNormal(model[i][0],model[i][1],model[i][2],normal);

        // change of coordinate
        
        for(int j=0;j<3;j++){
            
            float vector[4][1] = {model[i][j][0],model[i][j][1],model[i][j][2],1};

            float matrix[3][4] = {
                {eU[x],eU[y],eU[z],-dotProduct(OM,eU)},
                {eV[x],eV[y],eV[z],-dotProduct(OM,eV)},
                {eW[x],eW[y],eW[z],-dotProduct(OM,eW)},
            };

            float rslt[3][1];
            mulMat(3,4,matrix,4,1,vector,rslt);

            float matrix2[3][4] = {
                {Sx,0,0,0},
                {0,Sy,0,0},
                {0,0,(n+f)/(n-f),-1}      
            };
            
            
            float u_v_w_1[4][1] = {rslt[0][0],rslt[1][0],rslt[2][0],1} ;

            float rslt2[3][1];

            mulMat(3,4,matrix2,4,1,u_v_w_1,rslt2);
            
            
            for(int i_=0;i_<3;i_++){               
                model[i][j][i_]=((n-f)/(2*n*f*rslt[2][0]))*rslt2[i_][0];
                
            }

        }
        float cond = dotProduct(normal,L);
        
        if(cond < 0)shade[i] = -cond;
        else shade[i] = 0;
    }
    
    //allocate array W x H of pixel

    Pixel **picture = malloc((h+1000)*sizeof(Pixel*));
    assert(picture != NULL);
    for(int i=0;i<(h+1000);i++){
        picture[i]=malloc((w+1000)*sizeof(Pixel));
        assert(picture[i] != NULL);
    }
    //default value
    for(int i=0;i<(h+1000);i++){
        for (int j=0;j<(w+1000);j++){
            picture[i][j].R =CbgR;
            picture[i][j].G =CbgG;
            picture[i][j].B =CbgB;
            picture[i][j].depth = 1;
        }
    }
    
    //over each triangle
    startLoop = omp_get_wtime();
    #pragma omp parallel
    {   
        #pragma omp for schedule(dynamic)
        for(int t =0;t<T;t++){

            //over each Pixel
            float A[] = {model[t][a][x],model[t][a][y]};
            float B[] = {model[t][b][x],model[t][b][y]};
            float C[] = {model[t][c][x],model[t][c][y]};
            float x_min = A[x];
            float x_max = A[x];
            float y_min = A[y];
            float y_max = A[y];
            if(B[x] < x_min) x_min = B[x];
            if(C[x] < x_min) x_min = C[x];
            if(B[x] > x_max) x_max = B[x];
            if(C[x] > x_max) x_max = C[x];
            if(B[y] < y_min) y_min = B[y];
            if(C[y] < y_min) y_min = C[y];
            if(B[y] > y_max) y_max = B[y];
            if(C[y] > y_max) y_max = C[y];
            
            int i_start = floor(((float)h/2)*(1-y_max) -0.5);
            int i_end = ceil(((float)h/2)*(1-y_min) -0.5);
            int j_start = floor(((float)w/2)*(x_min+1)-0.5);
            int j_end = ceil(((float)w/2)*(x_max+1)-0.5);

            
            for(int i=i_start;i<i_end;i++){
                for (int j=j_start;j<j_end;j++){

                    //inside or not of the triangle
                    float P[] = { ( ( 2.*((float)j+0.5)) / (float)w ) - 1. , ( (-2.*((float)i+0.5)) / (float)h ) + 1. };
                    float A[] = {model[t][a][x],model[t][a][y]};
                    float B[] = {model[t][b][x],model[t][b][y]};
                    float C[] = {model[t][c][x],model[t][c][y]};
                    float AB[] = {B[x]-A[x],B[y]-A[y]};
                    float AP[] = {P[x]-A[x],P[y]-A[y]};
                    float AC[] = {C[x]-A[x],C[y]-A[y]};
                    float BC[] = {C[x]-B[x],C[y]-B[y]};
                    float BP[] = {P[x]-B[x],P[y]-B[y]};
                    float BA[] = {A[x]-B[x],A[y]-B[y]};

                    bool bool1 = ( AB[0]*AP[1] - AB[1]*AP[0] ) * ( AP[0]*AC[1] - AP[1]*AC[0] ) >= 0;
                    bool bool2 = ( BC[0]*BP[1] - BC[1]*BP[0] ) * ( BP[0]*BA[1] - BP[1]*BA[0] ) >= 0;

                    if (bool1 && bool2){ 
                        
                        float da = model[t][a][z];
                        float db = model[t][b][z];
                        float dc = model[t][c][z];

                        float wa = ((B[x]-A[x])*(P[y]-A[y])-(B[y]-A[y])*(P[x]-A[x]))/((B[x]-A[x])*(C[y]-A[y])-(B[y]-A[y])*(C[x]-A[x]));
                        float wb = ((A[x]-C[x])*(P[y]-C[y])-(A[y]-C[y])*(P[x]-C[x]))/((A[x]-C[x])*(B[y]-C[y])-(A[y]-C[y])*(B[x]-C[x]));
                        float wc = 1-wa-wb;

                        float dp = wa*da + wb*db + wc*dc;

                        if(dp < picture[i][j].depth){
                            
                            picture[i][j].R = shade[t]*CoR;
                            picture[i][j].G = shade[t]*CoG;
                            picture[i][j].B = shade[t]*CoB;
                            picture[i][j].depth = dp;
                        }
                    }
                }
            }
        }
    }
    endLoop = omp_get_wtime();

    //create and save the .ppm

    FILE *fp = fopen("image.ppm", "wb"); 
    
    fprintf(fp, "P6\n%d %d\n255\n",w,h);
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++){
            unsigned char color[3];
            color[0] = picture[i][j].R;  // red //
            color[1] = picture[i][j].G;  // green //
            color[2] = picture[i][j].B;  // blue //
            (void) fwrite(color,1,3,fp);
        }
    }
    fclose(fp);

    double end = omp_get_wtime(); //counter for total time of computation
    
    printf("time : %f\n",end - start);
    printf("time (loop): %f\n", endLoop - startLoop);

    for(int i =0;i<T;i++){
        for(int j=0;j<3;j++) free(model[i][j]); 
    }
    for(int i =0;i<T;i++) free(model[i]);
    free(model);
    free(shade);
    for(int i=0;i<(h+1000);i++)free(picture[i]);
    free(picture);
    
    return EXIT_SUCCESS;
}

float dotProduct(float vect_A[], float vect_B[]){
	float product = 0;
	// Loop for calculate dot product
	for (int i = 0; i < 3; i++)

		product = product + vect_A[i] * vect_B[i];
	return product;
}

void crossProduct(float vect_A[], float vect_B[], float cross_P[]){
	cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
	cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
	cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

void mulMat(int m_1,int n_1,float mat1[m_1][n_1],int m_2,int n_2,float mat2[m_2][n_2],float rslt[m_1][n_2]) {

    for (int i = 0; i < m_1; i++) {
        for (int j = 0; j < n_2; j++) {
            rslt[i][j] = 0;
            
            for (int k = 0; k < m_2; k++) {
                rslt[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

float norm(float vect[]){
    return sqrt(dotProduct(vect,vect));
}

void getNormal(float A[],float B[],float C[],float normal[]){
    
    float AB[] = {B[0]-A[0],B[1]-A[1],B[2]-A[2]};
    float AC[] = {C[0]-A[0],C[1]-A[1],C[2]-A[2]};
    float cross[3];
    crossProduct(AB,AC,cross);
    float norm_ = norm(cross);
    normal[0]=cross[0]/norm_;
    normal[1]=cross[1]/norm_;
    normal[2]=cross[2]/norm_;
}
