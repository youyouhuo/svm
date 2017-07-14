
#include "te.c"
//多项式核函数
double **getPolynomial(double **x,double **z,int row,int col,double c,int d){

double **x_transfer=getTransfer(x,row,col); //col*row

double **result=getTimes(x_transfer,z,col,row,col);


for(int i=0;i<col;++i){

	for(int j=0;j<col;++j){

		result[i][j]=pow(result[i][j]+c,d);
	}
}



freeMemory(x_transfer,col,row);


return result;
}

//求解核矩阵

double ** getKMatrix(double**(*func)(double **,double**,int,int,double,int ),double **input,double** t,int row,int col,double c,double d){

	double **result=getMemory(col,col);


	for(int i=0;i<col;++i){

		for(int j=0;j<col;++j){

//获取input中的某一列
			double **x_i=getOneCol(input,row,col,i);

			double **x_j=getOneCol(input,row,col,j);

//调用核函数计算K（x_i,x_j)
			double **calResult=func(x_i,x_j,row,1,c,d);


			freeMemory(x_i,row,1);

			freeMemory(x_j,row,1);

//计算(1/2)*t_i*t_j
			result[i][j]=0.5*t[i][0]*t[j][0]*calResult[0][0];


			freeMemory(calResult,1,1);


		}

	}


	return result;


}


//根据输入的核矩阵求解对应的导数矩阵

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

//对输入的线性方程组的矩阵进行LU分解
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


//根据LU分解得到的结果，求解线性方程组

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


//计算w*p

double **calculateWTimesP(double**(*func)(double **,double **,int,int,double,int),double **input,double **t,double **alpha,double **p,int row,int col,int pCol,int d,double c){



	double **result=getMemory(pCol,1);


	for(int i=0;i<pCol;++i){


		double **x_i=getOneCol(p,row,col,i);

		for(int index=0;index<col;++index){


			double **x_j=getOneCol(input,row,col,index);


			double **calResult=func(x_i,x_j,row,1,c,d);


			freeMemory(x_j,row,1);

			result[i][0]+=alpha[index][0]*t[index][0]*calResult[0][0];


			freeMemory(calResult,1,1);

		}


		freeMemory(x_i,row,1);



	}



	return result;






}


//计算得到偏置b
double getB(double**(*func)(double **,double **,int,int,double,int),double **input,double **t,double **alpha,double **p,int row,int col,int pCol,int d,double c){



	int none_zero_alpha=0;

	for(int i=0;i<col;++i){

		if(alpha[i][0]>0&&t[i][0]>0){

			none_zero_alpha=i;

			break;

		}
	}

	double b=0.0;

//选取alpha_i>0 并且 t_i>0 的正样本，计算b值

	if(none_zero_alpha<col&&t[none_zero_alpha][0]>0){

		double **x_i=getOneCol(input,row,col,none_zero_alpha);


		double **temp=calculateWTimesP(func,input,t,alpha,x_i,row,col,1,d,c);


		// printf(" none_zero_alpha:%d  alpha is :%lf temp is :%lf  t is %lf \n",none_zero_alpha,alpha[none_zero_alpha][0],temp[0][0],t[none_zero_alpha][0]);

		b=t[none_zero_alpha][0]-temp[0][0];


		freeMemory(x_i,row,1);

		freeMemory(temp,1,1);



	}


	return b;

}


//进行结果的预测

double **predict(double**(*func)(double **,double **,int,int,double,int),double **input,double **t,double **alpha,double **p,int row,int col,int pCol,int d,double c,double b){

//计算w*p

	double **result=calculateWTimesP(func,input,t,alpha,p,row,col,pCol,d,c);


	for(int i=0;i<pCol;++i){

		result[i][0]+=b;
	}


	return result;


}

//从文件中读取训练数据p，并且计算奇偶性，获取对应的t

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



void test1111(){


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


double c=1;
int d=2;

// int row=3;
// int col=8;

// double **input=getMemory(row,col);

// double **t=getMemory(col,1);

// getTrainningData("3.txt",input,col,row,t);   //p row*col



// // printArray(t,col,1);

// double c=1;
// int d=4;

// 	int row=6;
// 	int col=64;

// 	double **input=getMemory(row,col);

// 	double **t=getMemory(col,1);

// getTrainningData("6.txt",input,col,row,t);   //p row*col


// double c=1;
// int d=7;


printArray(t,col,1);




double ** k_matrix=getKMatrix(getPolynomial,input, t, row, col,c,d);  // col*col


// printf("111111");


double **dev_k_matrix=getDerivativeMatrix(k_matrix,col,col);  //col*(col+1)



// printf("1112111");

getLU(dev_k_matrix,col,col+1);



// printf("1113111");


double **alpha=getLUResult(dev_k_matrix,col,col+1);  //col*1


printf("input data  is:%dx%d\n",row,col);

printf("Polynomial:\n c is:%lf \n d is:%d\n",c,d);


printf("alpha is:\n");

printArray(alpha,col,1);


// printf("1114111");

double b=getB(getPolynomial,input,t,alpha,input,row,col,col,d, c);

double **result=predict(getPolynomial,input,t,alpha,input,row,col,col,d, c,b); //col*1



// printf("1115111");


// printArray(result,col,1);


getAddOrMinus(minus,result,t,col,1);
double totalError=0.0;

for(int i=0;i<col;++i)
	totalError+=fabs(result[i][0]);

printf("b is :%lf\n",b);


printf("total error is :%lf\n",totalError);




freeMemory(k_matrix,col,col);

freeMemory(dev_k_matrix,col,col+1);

freeMemory(alpha,col,1);

freeMemory(result,col,1);


freeMemory(t,col,1);
freeMemory(input,row,col);


}


int main(void){


	test1111();


	return 0;
}



