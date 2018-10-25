// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
// 
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice 
// Hall, John C. Hull,

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <cmath>


#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

// Multi-threaded pthreads header
#ifdef ENABLE_THREADS
#define MAX_THREADS 128
// Add the following line so that icc 9.0 is compatible with pthread lib.
#define __thread __threadp  
MAIN_ENV
#undef __thread
#endif

// Multi-threaded OpenMP header
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

// Multi-threaded header for Windows
#ifdef WIN32
#pragma warning(disable : 4305)
#pragma warning(disable : 4244)
#include <windows.h>
#define MAX_THREADS 128
#endif


//Precision to use for calculations
#define fptype float

#define NUM_RUNS 100
using namespace std;

typedef struct OptionData_ {
        fptype s;          // spot price
        fptype strike;     // strike price
        fptype r;          // risk-free interest rate
        fptype divq;       // dividend rate
        fptype v;          // volatility
        fptype t;          // time to maturity or option expiration in years 
                           //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)  
        char OptionType;   // Option type.  "P"=PUT, "C"=CALL
        fptype divs;       // dividend vals (not used in this test)
        fptype DGrefval;   // DerivaGem Reference Value
} OptionData;

OptionData *data;
fptype *prices;
int numOptions;

int    * otype;
fptype * sptprice;
fptype * strike;
fptype * rate;
fptype * volatility;
fptype * otime;
int numError = 0;
int nThreads;
fptype weight[6][6];
fptype error[1][6];

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

void getCofactor(fptype mat[6][6],fptype temp[6][6],int p, int q, int n)
{
	int i=0,j=0;
	for (int row =0; row<n; row++){
		for (int col =0; col<n; col++){
			if(row!=p && col!=q){
				temp[i][j++]=mat[row][col];

				if(j==n-1){
					j=0;
					i++;
				}
			}
		}
	}

}

fptype determinant(fptype mat[6][6],int n)
{
	fptype D=0;

	if(n==1)
		return mat[0][0];
	fptype temp[6][6];

	int sign =1;

	for(int f=0; f<n; f++)
	{
		getCofactor(mat,temp,0,f,n);
		D+=sign*mat[0][f]*determinant(temp,n-1);

		sign=-sign;
	}

	return D;
}

void adjoint(fptype A[6][6], fptype adj[6][6])
{
//	if(N==1)
//	{
//		adj[0][0]=1;
//		return;
//	}
	int sign =1;
	fptype temp[6][6];
	for(int i=1; i<6; i++)
	{
		for(int j=0; j<6; j++)
		{
			getCofactor(A,temp,i,j,6);

			sign =((i+j)%2==0)? 1:-1;

			adj[j][i]=(sign)*(determinant(temp,5));
		}
	}

}

bool inverse(fptype A[6][6],fptype inver[6][6])
{
	
	fptype det=determinant(A,6);
	if(det==0){
		cout<<"Singular matrix"<<endl;
		return false;
	}

	fptype adj[6][6];
	adjoint(A,adj);

	for(int i=0; i<6; i++)
		for(int j=0; j<6; j++)
			inver[i][j]=adj[i][j]/det;

	return true;

}

fptype CNDF ( fptype InputX ) 
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0) {
        InputX = -InputX;
        sign = 1;
    } else 
        sign = 0;

    xInput = InputX;
 
    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;
    
    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal   = xLocal_1 * xNPrimeofX;
    xLocal   = 1.0 - xLocal;

    OutputX  = xLocal;
    
    if (sign) {
        OutputX = 1.0 - OutputX;
    }
    
    return OutputX;
} 

// For debugging
void print_xmm(fptype in, char* s) {
    printf("%s: %f\n", s, in);
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/*
void error_injection(int *a,int rate) 
{
//	cout<<hex<<"original a="<<*a<<endl;
	unsigned int protect;//=0xfffff000;//e
	unsigned int error;//=  0x00000fff;//1
	switch (rate){
		case 0:{
			       protect=0xf0000000;
			       error  =0x0fffffff;
			       break;
		       }
		case 1:{
			       protect=0xff000000;
			       error  =0x00ffffff;
			       break;
		       }
		case 2:{
			       protect=0xfff00000;
			       error  =0x000fffff;
			       break;
		       }
		case 3:{
			       protect=0xffff0000;
			       error  =0x0000ffff;
			       break;
		       }
		case 4:{
			       protect=0xfffff000;
			       error  =0x00000fff;
			       break;
		       }
		case 5:{
			       protect=0xffffff00;
			       error  =0x000000ff;
			       break;
		       }
		case 6:{
			       protect=0xfffffff0;
			       error  =0x0000000f;
			       break;
		       }
		default:{//full protection
			       protect=0xffffffff;
			       error  =0x00000000;
			       break;
			}
	}
	float c=*a;
	int *d=(int *)&c;

	unsigned int upper=  0x00000000;//1
	*a&=protect;
	*d&=error;
	unsigned int num =upper;
	*d|=num;
	*a+=*d;
//	cout<<hex<<"approx a="<<*a<<endl<<endl;
}
*/

void error_injection(int *a,int rate) 
{
//	cout<<hex<<"original a="<<*a<<endl;
	unsigned int protect;//=0xfffff000;//e
	unsigned int error;//=  0x00000fff;//1
	switch (rate){
		case 0:{
			       protect=0xfff00000;
			       error  =0x000fffff;
			       break;
		       }
		case 1:{
			       protect=0xfffc0000;
			       error  =0x0003ffff;
			       break;
		       }
		case 2:{
			       protect=0xffff0000;
			       error  =0x0000ffff;
			       break;
		       }
		case 3:{
			       protect=0xffffc000;
			       error  =0x00003fff;
			       break;
		       }
		case 4:{
			       protect=0xfffff000;
			       error  =0x00000fff;
			       break;
		       }
		case 5:{
			       protect=0xfffffc00;
			       error  =0x000003ff;
			       break;
		       }
		case 6:{
			       protect=0xffffff00;
			       error  =0x000000ff;
			       break;
		       }
		case 7:{
			       protect=0xffffffc0;
			       error  =0x0000003f;
			       break;
		       }
		case 8:{
			       protect=0xfffffff0;
			       error  =0x0000000f;
			       break;
		       }
		case 9:{
			       protect=0xfffffffc;
			       error  =0x00000003;
			       break;
		       }
		
		default:{//full protection
			       protect=0xffffffff;
			       error  =0x00000000;
			       break;
			}
	}
	float c=*a;
	int *d=(int *)&c;

	unsigned int upper=  0x00000000;//1
	*a&=protect;
	*d&=error;
	unsigned int num =upper;
	*d|=num;
	*a+=*d;
//	cout<<hex<<"approx a="<<*a<<endl<<endl;
}



fptype approx(float a,int rate) 
{
	float c=a;
	error_injection((int *) &c,rate);
	return c;
}

fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
                            fptype strike, fptype rate, fptype volatility,
                            fptype time, int otype, float timet )
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1; 
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;  

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    xSqrtTime = sqrt(xTime);

    logValues = log( sptprice / strike );
        
    xLogTerm = logValues;
        
    
    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;
        
    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 -  xDen;

    d1 = xD1;
    d2 = xD2;
    
    NofXd1 = CNDF( d1 );
    NofXd2 = CNDF( d2 );

    FutureValueX = strike * ( exp( -(rate)*(time) ) );        
    if (otype == 0) {            
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else { 
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }
    
    return OptionPrice;
}
//**************************approximate function************************************//
fptype BlkSchlsEqEuroNoDiv_approx( fptype sptprice,
                            fptype strike, fptype rate, fptype volatility,
                            fptype time, int otype, float timet )
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1; 
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;    
   //error injection
//   error_injection((int *)&sptprice);
//   error_injection((int *)&strike);
//   error_injection((int *)&rate);
//final modification
//   error_injection((int *)&volatility);
//   error_injection((int *)&time);
//   error_injection((int *)&timet);

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    xSqrtTime = sqrt(xTime);

    logValues = log( sptprice / strike );
        
    xLogTerm = logValues;
        
    
    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;
        
    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 -  xDen;

    d1 = xD1;
    d2 = xD2;
    
    NofXd1 = CNDF( d1 );
    NofXd2 = CNDF( d2 );

    FutureValueX = strike * ( exp( -(rate)*(time) ) );        
    if (otype == 0) {            
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else { 
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }
    
    return OptionPrice;
}
/*
fptype error_predit(fptype in_error[5],fptype out_error,int iter, int protect[5]){

	fptype predit_out=1;
	if(iter<6){
		for(int i=0; i<iter;i++){
			int match=0;
			for (int j=0; j<5; j++){
				if(abs(in_error[j]==weight[iter][j])<0.1)
					match++;
			}
			cout<<"match ="<<match<<endl<<endl;
			if(match > 0){
//				iter--;
				return 0;
			}

		}
		error[0][iter]=out_error;
		for(int i=0; i<6; i++){
			if(i==5)
				weight[iter][i]=1;
			else
				weight[iter][i]=in_error[i];
		}
	}
	else{
		
		fptype inver[6][6];
		cout.precision(20);
		cout<<"weight = "<<endl;
		for(int i=0; i<6 ; i++){
			for(int j=0;j<6;j++){
				cout<<fixed<<weight[i][j]<<"  ";
			}
			cout<<endl;
		}
		cout<<"error ="<<endl;
		for(int i=0; i<6 ; i++){
			cout<<fixed<<error[0][i]<<endl;
		}
		cout<<"inverseable ="<<inverse(weight, inver)<<endl;
	
		
	}

	return predit_out;
}
*/
#define history_size 10
fptype error_threshold = 0.1;
fptype out_error_history[history_size];
int itemCount = 0;
int front = 0;
int rear = -1;
bool isFull(){
	return itemCount == history_size;
}
bool isEmpty(){
	return itemCount ==0;
}
void enqueue(fptype data){
	if(!isFull()){
		if(rear == history_size-1){
			rear = -1;
		}
		out_error_history[++rear] = data;
		itemCount++;
	}
}

fptype dequeue(){
	fptype data = out_error_history[front++];

	if(front == history_size){
		front = 0;
	}
	itemCount--;
	return data;
}


void error_config(fptype out_error, int protect[5]){
/*	for(int i=0;i<5;i++){
		cout<<"protect"<<i<<" = "<<protect[i]<<endl;
	}
	cout<<"----------------------"<<endl;
	cout<<"output_error = "<<out_error<<endl;
	cout<<"++++++++++++++++++++++"<<endl;
	cout<<endl;
*/
	if(!isFull()){
		enqueue(out_error);
	}
	else {
		dequeue();
		enqueue(out_error);
	}

	if(isFull()){
		fptype sum_his_error=0;
		for(int i=0; i<history_size; i++){
			sum_his_error=sum_his_error + out_error_history[i];
		}
		fptype average_his_error=sum_his_error / history_size;

		if(average_his_error >= error_threshold){ //exceed error threshold
			for(int i=0;i<5;i++){  //add protection
				protect[i]=protect[i]+2;
			}
		}
		else if((error_threshold - average_his_error) > 0.01){ //output error is below error threshold
			if(protect[0]>1){
				for(int i=0;i<5;i++){ //reduce protection
					protect[i]--;
				}
			}
		}
	}

}
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
#ifdef WIN32
DWORD WINAPI bs_thread(LPVOID tid_ptr){
#else
int bs_thread(void *tid_ptr) {
#endif
    int i, j, iter;
    fptype price;
    fptype approx_price;
    fptype priceDelta;
    int tid = *(int *)tid_ptr;
    int start = tid * (numOptions / nThreads);
    int end = start + (numOptions / nThreads);
    iter=0;
    int protect[5]={3,3,3,3,3};

    int re_execute=0;
    double sum_error=0;
//    int protect[][5];
    for (j=0; j<NUM_RUNS; j++) {
#ifdef ENABLE_OPENMP
#pragma omp parallel for
        for (i=0; i<numOptions; i++) {
#else  //ENABLE_OPENMP
        for (i=start; i<end; i++) {
#endif //ENABLE_OPENMP
            /* Calling main function to calculate option value based on 
             * Black & Sholes's equation.
             */
            price = BlkSchlsEqEuroNoDiv( sptprice[i], strike[i],
                                         rate[i], volatility[i], otime[i], 
                                         otype[i], 0);
//	    cout<<"sptprice ="<<sptprice[i]<<endl;
//	    cout<<"strike ="<<strike[i]<<endl;
//	    cout<<"rate ="<<rate[i]<<endl;
//	    cout<<"volatility ="<<volatility[i]<<endl;
//	    cout<<"otime ="<<otime[i]<<endl;
//	    cout<<"------------------------"<<endl;

	    fptype approx_sptprice=approx(sptprice[i],protect[0]);
	    fptype approx_strike=approx(strike[i],protect[1]);
	    fptype approx_rate=approx(rate[i],protect[2]);
	    fptype approx_volatility=approx(volatility[i],protect[3]);
	    fptype approx_otime=approx(otime[i],protect[4]);

//	    cout<<"approx sptprice ="<<approx_sptprice<<endl;
//	    cout<<"approx strike ="<<approx_strike<<endl;
//	    cout<<"approx rate ="<<approx_rate<<endl;
//	    cout<<"approx volatility ="<<approx_volatility<<endl;
//	    cout<<"approx otime ="<<approx_otime<<endl;
//	    cout<<"-------------------------"<<endl;

	    fptype error_sptprice=abs(sptprice[i]-approx_sptprice)/sptprice[i];
	    fptype error_strike=abs(strike[i]-approx_strike)/strike[i];
	    fptype error_rate=abs(rate[i]-approx_rate)/rate[i];
	    fptype error_volatility=abs(volatility[i]-approx_volatility)/volatility[i];
	    fptype error_otime=abs(otime[i]-approx_otime)/otime[i];

//	    cout<<"sptprice error="<<error_sptprice<<endl;
//	    cout<<"strike error="<<error_strike<<endl;
//	    cout<<"rate error="<<error_rate<<endl;
//	    cout<<"volatility error="<<error_volatility<<endl;
//	    cout<<"otime error="<<error_otime<<endl;

	    fptype in_error[5]= {error_sptprice,error_strike,error_rate,error_volatility,error_otime};

	    approx_price = BlkSchlsEqEuroNoDiv_approx( approx_sptprice, approx_strike,
                                         approx_rate, approx_volatility, approx_otime, 
                                         otype[i], 0);
	    fptype out_error;
	    if(price==0)
		    out_error=0;
	    else
		    out_error = abs( price - approx_price ) / price;
	    
	    sum_error = sum_error + out_error;

//	    cout<<"sum_error ="<<sum_error<<endl;
//	    cout<<"out_error ="<<out_error<<endl;
//	    cout<<endl;
//
//	    cout<<"price ="<<price<<endl;
//	    cout<<"approx_price ="<<approx_price<<endl;

	    error_config(out_error,protect);//config next cycle's approximation
	   
/*    	    if(out_error>=error_threshold){
		re_execute++;
		cout<<"reexecution times ="<<re_execute<<endl;
    	    }
*/
	    prices[i] = price;
            
//#ifdef ERR_CHK   
            priceDelta = data[i].DGrefval - price;
            if( fabs(priceDelta) >= 1e-4 ){
                printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
                       i, price, data[i].DGrefval, priceDelta);
                numError ++;
            }
//#endif
        }
	if(j==0){
//		cout<<"sum_error ="<<sum_error<<endl;
//		cout<<"end - start ="<<(end-start)<<endl;
		cout<<"total error ="<<sum_error/(end-start)<<endl;
		sum_error=0;
	}
    }

    return 0;
}

int main (int argc, char **argv)
{
    FILE *file;
    int i;
    int loopnum;
    fptype * buffer;
    int * buffer2;
    int rv;

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
        printf("PARSEC Benchmark Suite Version "__PARSEC_XSTRING(PARSEC_VERSION)"\n");
	fflush(NULL);
#else
        printf("PARSEC Benchmark Suite\n");
	fflush(NULL);
#endif //PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
   __parsec_bench_begin(__parsec_blackscholes);
#endif

   if (argc != 4)
        {
                printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
                exit(1);
        }
    nThreads = atoi(argv[1]);
    char *inputFile = argv[2];
    char *outputFile = argv[3];

    //Read input data from file
    file = fopen(inputFile, "r");
    if(file == NULL) {
      printf("ERROR: Unable to open file `%s'.\n", inputFile);
      exit(1);
    }
    rv = fscanf(file, "%i", &numOptions);
    if(rv != 1) {
      printf("ERROR: Unable to read from file `%s'.\n", inputFile);
      fclose(file);
      exit(1);
    }
    if(nThreads > numOptions) {
      printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");
      nThreads = numOptions;
    }

#if !defined(ENABLE_THREADS) && !defined(ENABLE_OPENMP)
    if(nThreads != 1) {
        printf("Error: <nthreads> must be 1 (serial version)\n");
        exit(1);
    }
#endif

    // alloc spaces for the option data
    data = (OptionData*)malloc(numOptions*sizeof(OptionData));
    prices = (fptype*)malloc(numOptions*sizeof(fptype));
    for ( loopnum = 0; loopnum < numOptions; ++ loopnum )
    {
        rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
        if(rv != 9) {
          printf("ERROR: Unable to read from file `%s'.\n", inputFile);
          fclose(file);
          exit(1);
        }
    }
    rv = fclose(file);
    if(rv != 0) {
      printf("ERROR: Unable to close file `%s'.\n", inputFile);
      exit(1);
    }

#ifdef ENABLE_THREADS
    MAIN_INITENV(,8000000,nThreads);
#endif
    printf("Num of Options: %d\n", numOptions);
    printf("Num of Runs: %d\n", NUM_RUNS);

#define PAD 256
#define LINESIZE 64

    buffer = (fptype *) malloc(5 * numOptions * sizeof(fptype) + PAD);
    sptprice = (fptype *) (((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = (int *) malloc(numOptions * sizeof(fptype) + PAD);
    otype = (int *) (((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

    for (i=0; i<numOptions; i++) {
        otype[i]      = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i]   = data[i].s;
        strike[i]     = data[i].strike;
        rate[i]       = data[i].r;
        volatility[i] = data[i].v;    
        otime[i]      = data[i].t;
    }

    printf("Size of data: %d\n", numOptions * (sizeof(OptionData) + sizeof(int)));

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_begin();
#endif
#ifdef ENABLE_THREADS
    int tids[nThreads];
    for(i=0; i<nThreads; i++) {
        tids[i]=i;
        CREATE_WITH_ARG(bs_thread, &tids[i]);
    }
    WAIT_FOR_END(nThreads);
#else//ENABLE_THREADS
#ifdef ENABLE_OPENMP
    {
        int tid=0;
        omp_set_num_threads(nThreads);
        bs_thread(&tid);
    }
#else //ENABLE_OPENMP
#ifdef WIN32 
    if (nThreads > 1)
    {
        HANDLE threads[MAX_THREADS];
                int nums[MAX_THREADS];
                for(i=0; i<nThreads; i++) {
                        nums[i] = i;
                        threads[i] = CreateThread(0, 0, bs_thread, &nums[i], 0, 0);
                }
                WaitForMultipleObjects(nThreads, threads, TRUE, INFINITE);
    } else
#endif
    {
        int tid=0;
        bs_thread(&tid);
    }
#endif //ENABLE_OPENMP
#endif //ENABLE_THREADS
#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_end();
#endif

    //Write prices to output file
    file = fopen(outputFile, "w");
    if(file == NULL) {
      printf("ERROR: Unable to open file `%s'.\n", outputFile);
      exit(1);
    }
    rv = fprintf(file, "%i\n", numOptions);
    if(rv < 0) {
      printf("ERROR: Unable to write to file `%s'.\n", outputFile);
      fclose(file);
      exit(1);
    }
    for(i=0; i<numOptions; i++) {
      rv = fprintf(file, "%.18f\n", prices[i]);
      if(rv < 0) {
        printf("ERROR: Unable to write to file `%s'.\n", outputFile);
        fclose(file);
        exit(1);
      }
    }
    rv = fclose(file);
    if(rv != 0) {
      printf("ERROR: Unable to close file `%s'.\n", outputFile);
      exit(1);
    }

#ifdef ERR_CHK
    printf("Num Errors: %d\n", numError);
#endif
    free(data);
    free(prices);

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_bench_end();
#endif

    return 0;
}

