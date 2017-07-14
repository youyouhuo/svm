#include "te.c"


double **getKernel_squre(double **input,int row,int col,double c,int *resultLength){


	int totalLength=1+row+row*row;

	*resultLength=totalLength;

	double **result=getMemory(totalLength,col);



	for(int colIndex=0;colIndex<col;++colIndex){


		result[0][colIndex]=c;

		int resultRow=1;

		double countNum=sqrt(2*c);

		for(int rowIndex=0;rowIndex<row;++rowIndex){

			result[resultRow][colIndex]=countNum*input[rowIndex][colIndex];

			++resultRow;
		}


		for(int rowIndex=0;rowIndex<row;++rowIndex){


			for(int rowIndex2=0;rowIndex2<row;++rowIndex2){

				result[resultRow][colIndex]=input[rowIndex][colIndex]*input[rowIndex2][colIndex];
				++resultRow;

			}

		}


	}

	return result;
}





double ** getKMatrix(double**(*func)(double **,int,int,double,int *),double **input,double** t,int row,int col){


	// int totalLength=1+row+row*row;

	double **result=getMemory(col,col);


	for(int i=0;i<col;++i){



		for(int j=0;j<col;++j){


			double **x_i=getOneCol(input,row,col,i);

			double **x_j=getOneCol(input,row,col,j);

			int newRowLength=0;

			double **theta_i=func(x_i,row,1,1,&newRowLength);

			newRowLength=0;

			double **theta_j=func(x_j,row,1,1,&newRowLength);


			freeMemory(x_i,row,1);

			freeMemory(x_j,row,1);


			for(int calIndex=0;calIndex<newRowLength;++calIndex){

				result[i][j]=result[i][j]+theta_i[calIndex][0]*theta_j[calIndex][0];


			}
// printf("result :%lf \n",result[i][j]);

			result[i][j]=0.5*t[i][0]*t[j][0]*result[i][j];

			// printf("result :%lf \n",result[i][j]);


			freeMemory(theta_i,newRowLength,1);

			freeMemory(theta_j,newRowLength,1);



		}
	
	}


	return result;


}


double ** getDerivativeMatrix(double **input,int row,int col){

double **result=getMemory(row,col+1);


for(int i=0;i<row;++i){

	result[i][col]=1;
}

for(int i=0;i<row;++i){

	for(int j=0;j<col;++j){

		result[i][j]=input[i][j]+input[j][i];
	}
}


return result;
}


void getLU(double **input,int row,int col){


for(int index=0;index<row;++index){


for(int tempCol=index;tempCol<col;++tempCol){

	for(int i=0;i<index;++i){
input[index][tempCol]-=input[index][i]*input[i][tempCol];

	}
	
}


for(int tempRow=index+1;tempRow<row;++tempRow){

for(int i=0;i<index;++i){

	input[tempRow][index]-=input[tempRow][i]*input[i][index];
}


input[tempRow][index]=input[tempRow][index]/input[index][index];

}

}

}


double **getLUResult(double **input,int row,int col){

double **result=getMemory(row,1);

for(int i=row-1;i>=0;--i){

double temp=input[i][col-1];

for(int k=i+1;k<row;++k){

	temp=temp-result[k][0]*input[i][k];
}

result[i][0]=temp/input[i][i];



}


return result;

}



double **getWeight(double**(*func)(double **,int,int,double,int *),double **input,double **alpha,double **t,int row,int col,int *Length,double *b){

double **weight=NULL;


for(int i=0;i<col;++i){

double **x_i=getOneCol(input,row,col,i);

int totalLength=0;
double **theta_i=func(x_i,row,1,1,&totalLength);


*Length=totalLength;

// printArray(theta_i,totalLength,1);


freeMemory(x_i,row,1);



if(weight==NULL){

	weight=getMemory(totalLength,1);
}

for(int index=0;index<totalLength;++index){

	weight[index][0]+=t[i][0]*alpha[i][0]*theta_i[index][0];
}

freeMemory(theta_i,totalLength,1);


}


int none_zero_alpha=0;

for(int i=0;i<col;++i){

	if(alpha[i][0]>0){
    
    none_zero_alpha=i;

    break;
	
	}
}


if(none_zero_alpha<col){


double **x_i=getOneCol(input,row,col,none_zero_alpha);

int totalLength=0;

double **theta_i=func(x_i,row,1,1,&totalLength);

double temp=0.0;

for(int i=0;i<totalLength;++i){

	temp+=weight[i][0]*theta_i[i][0];
}




*b=t[none_zero_alpha][0]-temp;


// printf("t is:%lf  b is:%lf temp is :%lf  alpha is :%lf\n",t[none_zero_alpha][0],*b,temp,alpha[none_zero_alpha][0]);

freeMemory(x_i,row,1);

freeMemory(theta_i,totalLength,1);


}





return weight;

}



void getBitArray(double **array,int n,int number){

	if(number>=n){


printArray(array,1,n);
return;

	}


	array[0][number]=-1;
		
	getBitArray(array,n,number+1);

	array[0][number]=1;

getBitArray(array,n,number+1);


}



void testInput(){

double **input=getMemory(2,4);

int row=2;
int col=4;

input[0][0]=-1;
input[1][0]=-1;

input[0][1]=-1;
input[1][1]=1;

input[0][2]=1;
input[1][2]=-1;

input[0][3]=1;
input[1][3]=1;


double **t=getMemory(4,1);

t[0][0]=-1;

t[1][0]=1;
t[2][0]=1;
t[3][0]=-1;


double **result=getKMatrix(getKernel_squre,input,t,2,4);

// printArray(result,4,4);

double **result1=getDerivativeMatrix(result,4,4);
// 
// printArray(result1,4,5);

getLU(result1,4,5);

double **me=getLUResult(result1,4,5);

// printArray(me,4,1);



int totalLength=0;

double b=0.0;

double **weight=getWeight(getKernel_squre,input,me,t,2,4, &totalLength, &b);




printArray(weight,totalLength,1);






double **calInput=getKernel_squre(input,2,4,1,&totalLength);

double **weightTranfer=getTransfer(weight,totalLength,1);


double **calresult=getTimes(weightTranfer,calInput,1,totalLength,4);


printf("result is:\n");

for(int ou=0;ou<4;++ou){

printf("%lf ",calresult[0][ou]);

}




printf("\n");







// double b=0.0;

// double cal=0.0;

// double **the

// for(int i=0;i<totalLength;++i){

// cal+=weight[i][0]*


// }



freeMemory(weight,totalLength,1);

freeMemory(me,4,1);

freeMemory(result,4,4);

freeMemory(result1,4,5);


freeMemory(t,4,1);
freeMemory(input,2,4);



}



void testLU(){



int row=4;int col=5;
double **temp=getMemory(row,col);

temp[0][0]=1;
temp[0][1]=2;
temp[0][2]=-12;
temp[0][3]=8;
temp[0][4]=27;

temp[1][0]=5;
temp[1][1]=4;
temp[1][2]=7;
temp[1][3]=-2;
temp[1][4]=4;


temp[2][0]=-3;
temp[2][1]=7;
temp[2][2]=9;
temp[2][3]=5;
temp[2][4]=11;

temp[3][0]=6;
temp[3][1]=-12;
temp[3][2]=-8;
temp[3][3]=3;
temp[3][4]=49;


// printArray(temp,3,4);

getLU(temp,row,col);

// printArray(temp,3,4);



double **re=getLUResult(temp,row,col);

printArray(re,row,1);

freeMemory(re,row,1);
freeMemory(temp,row,col);





}








void testKernel(){

double **temp=getMemory(4,4);

double num=1;

for(int i=0;i<4;++i){


	for(int j=0;j<4;++j){

temp[j][i]=num;
num=num+1;

	}
}


printArray(temp,4,4);

int resultLength=0;

double **result=getKernel_squre(temp,4,4,1,&resultLength);

printArray(result,21,4);


printf("result Length:%d\n",resultLength);

freeMemory(temp,4,4);


freeMemory(result,21,4);




}


void testBitArray(){

int row=1;
int col=6;

double **array=getMemory(row,col);

// for(int i=0;i<col;++i){

// 	array[0][i]=-1;
// }


getBitArray(array,col,0);

freeMemory(array,row,col);



}



void getTrainningData(char*path,double **p,int row,int col,double **t){


	double **p_data=loadData(path,row,col);

	int colIndex=0;

	for(int i=0;i<row;++i){

		int count=0;

		for(int j=0;j<col;++j){

           p[j][i]=p_data[i][j];

           if(p_data[i][j]>0) ++count;


		}


		if(count%2==0){
			t[colIndex][0]=-1;
		}else{

			t[colIndex][0]=1;
		}

		++colIndex;
	

	}




}


void testLoadFile(){

int row=6;
int col=64;

double **p=getMemory(row,col);

double **t=getMemory(col,1);

getTrainningData("6.txt",p,col,row,t);

printArray(p,row,col);


printArray(t,col,1);


}


void test_of_two(){


int row=2;
int col=4;


double **p=getMemory(row,col);

double **t=getMemory(col,1);


p[0][0]=-1;
p[1][0]=-1;

p[0][1]=-1;
p[1][1]=1;

p[0][2]=1;
p[1][2]=-1;

p[0][3]=1;
p[1][3]=1;


t[0][0]=-1;

t[1][0]=1;
t[2][0]=1;
t[3][0]=-1;


double **result=getKMatrix(getKernel_squre,p,t,row,col);  // col*col


double **result1=getDerivativeMatrix(result,col,col);   //col*(col+1)


freeMemory(result,col,col);



getLU(result1,col,col+1);

// printArray(result1,col,col+1);

double **lu_result=getLUResult(result1,col,col+1);  //col*1


freeMemory(result1,col,col+1);


int totalLength=0;

double b=0.0;

double **weight=getWeight(getKernel_squre,p,lu_result,t,row,col, &totalLength, &b);  //totalLength*1


freeMemory(lu_result,col,1);



double **calInput=getKernel_squre(p,row,col,1,&totalLength);

double **weightTranfer=getTransfer(weight,totalLength,1);

double **calresult=getTimes(weightTranfer,calInput,1,totalLength,col);


printf("result is:\n");

for(int ou=0;ou<col;++ou){

printf("%lf ",calresult[0][ou]+b);

}




printf("\n");


freeMemory(weight,totalLength,1);

freeMemory(weightTranfer,1,totalLength);

freeMemory(calInput,totalLength,col);

freeMemory(calresult,1,col);

freeMemory(p,row,col);

freeMemory(t,col,1);


}



void test_of_three(){

int row=3;
int col=8;

double **p=getMemory(row,col);

double **t=getMemory(col,1);

getTrainningData("3.txt",p,col,row,t);   //p row*col


// printArray(p,row,col);

// printArray(t,col,1);


double **result=getKMatrix(getKernel_squre,p,t,row,col);  // col*col

// printArray(result,col,col);

double **result1=getDerivativeMatrix(result,col,col);   //col*(col+1)


// printArray(result1,col,col+1);


freeMemory(result,col,col);




getLU(result1,col,col+1);

printArray(result1,col,col+1);

double **lu_result=getLUResult(result1,col,col+1);  //col*1


freeMemory(result1,col,col+1);


// printArray(lu_result,col,1);

int totalLength=0;

double b=0.0;

double **weight=getWeight(getKernel_squre,p,lu_result,t,row,col, &totalLength, &b);  //totalLength*1


freeMemory(lu_result,col,1);



double **calInput=getKernel_squre(p,row,col,1,&totalLength);

double **weightTranfer=getTransfer(weight,totalLength,1);

double **calresult=getTimes(weightTranfer,calInput,1,totalLength,col);


printf("result is:\n");

for(int ou=0;ou<col;++ou){

printf("%lf ",calresult[0][ou]+b);

}




printf("\n");


freeMemory(weight,totalLength,1);

freeMemory(weightTranfer,1,totalLength);

freeMemory(calInput,totalLength,col);

freeMemory(calresult,1,col);

freeMemory(p,row,col);

freeMemory(t,col,1);


}




int main(void){

// testKernel();

// testLU();

	// testInput();

// testBitArray();


	// testLoadFile();


	// test_of_three();

	test_of_two();


	return 0;
}



