	#include<stdio.h>
	#include<malloc.h>
#include <math.h>

#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double **getMemory(int row,int col){
	
	double **result=(double**)malloc(row*sizeof(double*));
	
	int i=0;
	for(i=0;i<row;++i){

		result[i]=(double*)malloc(col*sizeof(double));
	}
	
	int j=0;
	
	for(i=0;i<row;++i){
		for(j=0;j<col;++j)
		{

			result[i][j]=0.0;


		}
	}
	
	return result;
}

void freeMemory(double **input,int row,int col){

	for(int i=0;i<row;++i){



		free(input[i]);
	}

	free(input);

}


//产生随机数
int getRandom(int max){

 int number=rand()%max;

return number;

}
//输出矩阵
void printArray(double **outData,int row,int col){
	
	int i=0;
	int j=0;
	for( i=0;i<row;++i){

		for(j=0;j<col;j++){

			printf("%lf ",outData[i][j]);
		}

		printf("\n");
	}
	
}

//矩阵乘法

double** getTimes(double **leftOne,double **right,int leftOneRow,int leftOneCol,int rightOneCol){
	
	int i=0;
	int j=0;
	int k=0;
	double ** result=getMemory(leftOneRow,rightOneCol);
	
	for(i=0;i<leftOneRow;i++){

		for(j=0;j<rightOneCol;j++){

			for(k=0;k<leftOneCol;k++){

				result[i][j]+=leftOne[i][k]*right[k][j];

			}
		}
	}
	
	return result;
}


//矩阵转置
double **getTransfer(double **input,int row,int col){
	
double **result=getMemory(col,row);
	for(int i=0;i<row;++i){

		for(int j=0;j<col;++j){

			result[j][i]=input[i][j];
		}
	}

	return result;
}

//代数余子式

double **getCofactor(double **input,int row,int col,int rowNumber,int colNumber){


if(row<=0) return NULL;
double **result=getMemory(row-1,col-1);

int i0=0;
int j0=0;

for(int i=0;i<row;++i){

if(i==rowNumber) continue;

j0=0;
	for(int j=0;j<col;++j){

if(j==colNumber) continue;

result[i0][j0]=input[i][j];

j0++;
	}
	i0++;
}

return result;


}

// 行列式

double getDeterminant(double **input,int row,int col){


if(row==2) return input[0][0]*input[1][1]-input[0][1]*input[1][0];

if(row==1) return input[0][0];




double sum=0.0;

for(int i=0;i<col;++i){

double **tempInput=getCofactor(input,row,col,0,i);
double tempSum=input[0][i]*getDeterminant(tempInput,row-1,col-1);

	if(i%2==0) sum+=tempSum;
	else sum+=-tempSum;

	freeMemory(tempInput,row-1,col-1);
}

return sum;

}

//矩阵的逆矩阵
double **getINverse(double **input,int row,int col){

double **result=getMemory(row,row);

double Determinant=getDeterminant(input,row,row);

for(int i=0;i<row;++i){

	for(int j=0;j<row;++j){

		double **tempCofactor=getCofactor(input,row,row,i,j);

    if((i+j)%2==0){

    	result[j][i]=getDeterminant(tempCofactor,row-1,row-1)/Determinant;

    }else result[j][i]=-getDeterminant(tempCofactor,row-1,row-1)/Determinant;


freeMemory(tempCofactor,row-1,row-1);

	}
}



return result;

}

//两个矩阵的对位操作

void getAddOrMinus(double (*fun)(double,double),double** firstOne,double **secondOne,int row,int col){

for(int i=0;i<row;++i){

for(int j=0;j<col;++j){

firstOne[i][j]=fun(firstOne[i][j],secondOne[i][j]);

}

}

}

//硬极限函数
double hardlim(double input){


	if(input<0.0) return -1;
	else return 1;
}


//两个数相加
double add(double first,double second){

	return first+second;
}
//两个数相减
double minus(double first,double second){

return first-second;

}
  
//将对矩阵中的每一个数调用fun函数
void activateMatrix(double(*fun)(double),double** input,int row,int col){

for(int i=0;i<row;++i){

	for(int j=0;j<col;++j){

		input[i][j]=fun(input[i][j]);

	}
}

}


//从指定文件路径中加载数据
double **loadData(char *filePath,int row,int col){

if(filePath==NULL) return NULL;

double **result=getMemory(row, col);

FILE *fp=fopen(filePath,"r");

if(fp==NULL) return result; 

for(int i=0;i<row;++i){
	for(int j=0;j<col;++j){
		fscanf(fp,"%lf",&result[i][j]);
	}
}

fclose(fp);

return result;

}




//数据矩阵对应的图形
void printChart(double **input,int row,int col){


	for(int i=0;i<row;++i){

if(i%5==0) printf("\n");

			if(input[i][0]>0.0)
			printf("*");
		else printf(" ");

		
	}


	printf("\n");
}

//将数据保存到文件中
void saveData(char *filePath,double **result,int row,int col){

if(filePath==NULL) return;

FILE *fp=fopen(filePath,"w");

if(fp==NULL) return ; 

for(int i=0;i<row;++i){
	for(int j=0;j<col;++j){
		fprintf(fp,"%lf ",result[i][j]);
	}
	fprintf(fp, "\n");
}

fclose(fp);


}

//获取矩阵中的某一行
double **getOneRow(double **input,int row,int col,int rowNumber){

double **result=getMemory(col,1);

for(int i=0;i<col;++i){

	result[i][0]=input[rowNumber][i];
}


return result;

}
//获取矩阵中的某一列

double **getOneCol(double **input,int row,int col,int colNumber){

double **result=getMemory(row,1);

for(int i=0;i<row;++i){

	result[i][0]=input[i][colNumber];
}


return result;

}

//训练数据  w=t*(p_transfer*p)^-1 *p_transfer
double **trainningData(double **p_transfer,double **output,int inputRow,int inputCol,int outputRow,int outputCol ){


double **p=getTransfer(p_transfer,inputRow,inputCol);  //30*10

double **p_transfer_times_p=getTimes(p_transfer,p,inputRow,inputCol,inputRow);  //10*30*30*10=10*10

double **inv=getINverse(p_transfer_times_p,inputRow,inputRow);  //10*10

double **p_plus=getTimes(inv,p_transfer,inputRow,inputRow,inputCol);  //10*10 *10*30=10*30

double **weight=getTimes(output,p_plus,outputRow,outputCol,inputCol); //30*10 *10*30=30*30

// saveData("weight",weight,outputRow,inputCol);


freeMemory(p,inputCol,inputRow);

freeMemory(p_transfer_times_p,inputRow,inputRow);


freeMemory(inv,inputRow,inputRow);

freeMemory(p_plus,inputRow,inputCol);

return weight;






}

void printResult(double **weight,double **result,double **p,double **t_transfer){

double **result_transfer=getTransfer(result,30,10);


double error1=0.0;

for(int i=0;i<10;++i){

	double **oneTem=getOneRow(result_transfer,10,30,i);

	double **oneInput=getOneCol(p,30,10,i);

	double tempError=0.0;

	for(int j=0;j<30;++j){

		tempError+=fabs(oneTem[j][0]-t_transfer[i][j]);
	}

	error1+=tempError;

	printf("输入数据如下：\n");

	printChart(oneInput,30,1);

	printf("结果数据如下：\n");

	printChart(oneTem,30,1);

	freeMemory(oneTem,30,1);

	freeMemory(oneInput,30,1);

	printf("tempError is %lf \n",tempError);


}


printf("total error is :%lf \n",error1);


freeMemory(result_transfer,10,30);


}


void testTrain(){

double **p_transfer=loadData("da.txt",10,30);

double **t_transfer=loadData("re.txt",10,30);

double **t=getTransfer(t_transfer,10,30);



double **weight=trainningData(p_transfer,t,10,30,30,10);



double **p=getTransfer(p_transfer,10,30);

double **result=getTimes(weight,p,30,30,10);

activateMatrix(hardlim,result,30,10);

printResult(weight,result,p,t_transfer);


freeMemory(result,30,10);

result=NULL;


printf("7位随机噪声输入时：\n");

for(int i=0;i<10;++i){

for(int j=0;j<7;++j){

int index=getRandom(30);

	if(p[index][i]>0.0) p[index][i]=-1.0;

else p[index][i]=1.0;
}

}




result=getTimes(weight,p,30,30,10);

activateMatrix(hardlim,result,30,10);

printResult(weight,result,p,t_transfer);






freeMemory(result,30,10);

freeMemory(t,30,10);

freeMemory(t_transfer,10,30);

freeMemory(p_transfer,10,30);

freeMemory(p,30,10);



freeMemory(weight,30,30);


}




void testMHLS(){

double **input=loadData("da.txt",10,30);

double **input_transfer=getTransfer(input,10,30);

double **times=getTimes(input,input_transfer,10,30,10);

double result=getDeterminant(times,10,10);

printf("------%lf====\n",result);


double **inv=getINverse(times,10,10);


saveData("com",inv,10,10);

printArray(inv,10,10);

freeMemory(times,10,10);

freeMemory(inv,10,10);
freeMemory(input_transfer,30,10);




freeMemory(input,10,30);

}



void testSelf(){
double **input=loadData("sef.txt",9,30);

double **p_transfer=getMemory(3,30);




for(int i=0;i<3;++i){

	for(int j=0;j<30;++j){

		p_transfer[i][j]=input[i][j];
	}


}

double **t=getTransfer(p_transfer,3,30);

double **weight=getTimes(t,p_transfer,30,3,30);

char *inferData[4];

inferData[0]="100%输入时：";

inferData[1]="50%输入时：";


inferData[2]="33%输入时：";

inferData[3]="7位随机噪声输入时：";



for(int i=0;i<9;++i){

double**oneTestData=getOneRow(input,9,30,i);

printf("%s\n",inferData[(i)/3]);

printf("输入数据如下：\n");

printChart(oneTestData,30,1);

double **result=getTimes(weight,oneTestData,30,30,1);

activateMatrix(hardlim,result,30,1);

printf("结果数据如下：\n");

printChart(result,30,1);

printf("\n");


getchar();

freeMemory(result,30,1);

freeMemory(oneTestData,30,1);

}



for(int i=0;i<3;++i){

for(int j=0;j<7;++j){

int index=getRandom(30);

// printf("%d \n",index);

if(p_transfer[i][index]>0.0) p_transfer[i][index]=-1.0;

else p_transfer[i][index]=1.0;


}

}

for(int i=0;i<3;++i){

double**oneTestData=getOneRow(p_transfer,3,30,i);

printf("%s\n",inferData[3]);

printf("输入数据如下：\n");

printChart(oneTestData,30,1);

double **result=getTimes(weight,oneTestData,30,30,1);

activateMatrix(hardlim,result,30,1);

printf("结果数据如下：\n");

printChart(result,30,1);

printf("\n");


getchar();
freeMemory(result,30,1);

freeMemory(oneTestData,30,1);

}











freeMemory(p_transfer,3,30);

freeMemory(t,30,3);


freeMemory(input,9,30);

}



void test(){


int row=4;
int col=6;
	double **res=getMemory(row,col);

	double **res2=getMemory(row,col);

	int i=0;

	int j=0;

	

	double temp=1.0;
	for( i=0;i<row;++i){

		for(j=0;j<col;j++){

			res[i][j]=temp;
			res2[i][j]=temp;

			temp=temp+1;

		}

	}
	
	// printArray(res,4,4);
	// printArray(res2,4,4);
	
	
	// double **resNow=getTimes(res,res2,4,4,4);
	
	printArray(res,row,col);

	double **trander=getTransfer(res,row,col);

	// double **resNow=getTimes(res,trander,row,col,row);

	double **resNow;
    resNow=getTimes(res,trander,row,col,row);

	printArray(trander,col,row);

	printArray(resNow,row,row);

	freeMemory(trander,col,row);

	freeMemory(resNow,row,row);

	freeMemory(res,row,col);
	freeMemory(res2,row,col);

}



int main1(void){

	srand( (unsigned)time( NULL ) );

	// testSelf();

	testTrain();
	

	return 0;
} 
