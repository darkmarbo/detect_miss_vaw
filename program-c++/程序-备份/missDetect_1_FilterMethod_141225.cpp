
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>



#define PI (3.1415926535897932384626433832795)

#define posThres 600 //13.6ms under 44.1kHz
#define valThres 10
#define varThres 0.005

#define SAMPLING_FREQENCY 44100
#define MAX_SIG_LENGTH (SAMPLING_FREQENCY*10)

float Input[MAX_SIG_LENGTH];
float output[MAX_SIG_LENGTH];
long sigLength;

typedef struct _wavhead
{
	char            riff[4];            //"RIFF"
	unsigned long   filelong;           // +8 = File size
	char            wav[8];             //"WAVEfmt "
	unsigned long   t1;                 
	short           tag;
	short           channels;
	unsigned long   samplerate;         
	unsigned long   typesps;            
	unsigned short  psbytes;            
	unsigned short  psbits;             
	char            data[4];            
	unsigned long   sumbytes;           
}WAVEHEAD;

short *allbuf ;//全部数据

int ReadFile(char *wfile)
{
	bool oflag=false;
	FILE *fp=NULL;
	WAVEHEAD head;
	int SAMFREQ=-1;
	int sample_count=0,channel_num=0,readflag=0;
	try
	{
		//判断声音文件
		if (strstr(wfile, ".wav")) {
			fp=fopen(wfile, "rb");
			if (fp == NULL) {
				return -1;
			}
			oflag=true;
			fseek(fp,0,SEEK_END);
			sample_count = ftell(fp) - sizeof(WAVEHEAD);
			fseek(fp,0,SEEK_SET);
			fread(&head, 1, sizeof(WAVEHEAD), fp);
			//data
			if(head.data[0]!='d'&&head.data[1]!='a'&&head.data[2]!='t'&&head.data[3]!='a')
			{
				fclose(fp);
				return -1;
			}
			//RIFF
			if(head.riff[0]!='R'&&head.riff[1]!='I'&&head.riff[2]!='F'&&head.riff[3]!='F')
			{
				fclose(fp);
				return -1;
			}
			//"WAVEfmt "
			if(head.wav[0]!='W'&&head.wav[1]!='A'&&head.wav[2]!='V'&&head.wav[3]!='E'&&head.wav[4]!='f'&&head.wav[5]!='m'&&head.wav[6]!='t'&&head.wav[7]!=' ')
			{
				fclose(fp);
				return -1;
			}
			//定位数据
			fseek(fp,(long)(head.t1-16)-4,SEEK_CUR);
			fread(&head.sumbytes,1,sizeof(long),fp);
			//得到字节数
			sample_count=head.sumbytes;
			if(head.samplerate>48000||head.samplerate<0)
			{
				fclose(fp);
				exit(-1);
			}
			SAMFREQ = head.samplerate;
			channel_num = head.channels;
		}
		//得到样本数（n个通道样本数和，且为16bit）
		sample_count /= sizeof(short);
		if (sample_count % channel_num != 0) {
			fclose(fp);
			return -2;
		}
		//分配空间读取数据
		allbuf = (short*)malloc(sample_count * sizeof(short));
		fread(allbuf, sizeof(short), sample_count,fp);
		fclose(fp);
		oflag=false;

		return sample_count;
	}
	catch(...)
	{
		if(oflag)
			fclose(fp);

		if(allbuf)free(allbuf);
		allbuf=NULL;
		return -1;

	}
}



float bessi0(float x)
{
	float ax,ans;
	float y;

	ax = fabs(x);
	if (ax < 3.75)
	{
		y = x/3.75;
		y *= y;
		ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y*0.45813e-2)))));
	}
	else
	{
		y = 3.75 / ax;
		ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
	}
	return ans;
}

void kaiserPara(float delta, float transBw, int &fLength, float &beta)
{
	float a, len;
	
	a= -20 * log10(delta);

	if (a <= 21) beta = 0;
	else if (a<= 50) beta = 0.5842 * pow(float(a-21), float(0.4)) + 0.07889 * (a-21);
	else beta = 0.1102 * (a - 8.7);

	len = (a - 7.95) / 14.36 / transBw;
	fLength = int(len);
	if ((len - fLength) < 0.5) fLength++;
	else fLength+=2;

	if (fLength%2 != 0) fLength++;
}

void kaiserLowPass(float *filter, int fLength, float beta, float wn)
{
	int tim, step;
	float k, sum;

	for (tim=0; tim<=fLength; tim++)
	{
		k = 2*tim/float(fLength) - 1;
		filter[tim] = bessi0( beta*sqrt(1- k*k)) / bessi0( beta );
	}

	sum=0;
	for (tim=0; tim<=fLength; tim++)
	{
		step = tim - fLength/2;
		if (step !=0) filter[tim] *= sin(wn * PI * step) / PI / step;
		else filter[tim] *= wn;

		sum += filter[tim];
	}
}

void kaiserBandPass(float *filter, int fLength, float beta, float wc, float wn)
{
	int tim;

	kaiserLowPass(filter, fLength, beta, wn);

	for (tim=0; tim<=fLength; tim++)
		filter[tim] *= 2 * cos((tim-fLength/2) * wc * PI);
}
float * leapDetect(char * filename,int &len);
void main(int argc, char *argv[])
{
	char* filename = argv[1];
	int len;
    float* leapPositions = leapDetect(filename,len);
	if (leapPositions==NULL)
	{
		printf("read error...");
		return;
	}

	FILE* fp;
	if((fp=fopen(argv[2],"a+"))==NULL)
	{
		printf("Please provide the second parameter...");
		return;
	}
	else
	{
		fprintf(fp,"%d possible jumping positions:\n",len);
		for (int i=0;i<len;i++)
		{			
			fprintf(fp,"%.2fs ",leapPositions[i]);
		}
		fprintf(fp,"\n");
	}
	
	/*printf("There are %d possible jumping positions laocated at:\n",len);
	for (int i=0;i<len;i++)
	{
		printf("%.2fs ",leapPositions[i]);
	}
	printf("\n");
	*/
	

	
}

float * leapDetect(char * filename,int &len)
{
	float *filter;
	float beta;
	int fLength;

	sigLength = ReadFile(filename);
	if (sigLength<1)
	{
		return NULL;
	}
	for (int i=0;i<sigLength;i++)
	{
		Input[i] = allbuf[i]/float(32768);
	}

	kaiserPara(0.01, 1000/float(SAMPLING_FREQENCY), fLength, beta);		
	filter = new float[fLength+1];
	kaiserBandPass(filter, fLength, beta, 40050/float(SAMPLING_FREQENCY), 4050/float(SAMPLING_FREQENCY));

	for(long n=0; n<sigLength; n++)
	{
		for(int m=0; m<=fLength; m++)
		{				
			long tim = n + fLength/2 - m;
			if ( (tim >= 0) && (tim < sigLength) ) 
				output[n] += Input[tim] * filter[m];
		}
	}
	delete filter;

	//mean
	float mu = 0;
	for (int i=0;i<sigLength; i++)
	{
		mu += output[i];
	}
	mu /= float(sigLength);
	//std
	float std = 0;
	for (int i=0;i<sigLength;i++)
	{
		std += pow(float(output[i]-mu),float(2.0));
	}
	std = pow(std/float(sigLength),float(0.5));
	//normalization
	float* output_norm = (float*)malloc(sigLength * sizeof(float));
	for (int i=0;i<sigLength;i++)
	{
		output_norm[i] = (output[i]-mu)/float(std);
	}
	float* output_norm_select = (float*)malloc(sigLength * sizeof(float));
	//pick out candidate points
	for (int i=0;i<sigLength;i++)
	{
		if (73620==i)
		{
			printf("");
		}
		if (output_norm[i]> valThres && output_norm[i]*std > varThres )
		{
			output_norm_select[i] = output_norm[i];
		}
		else
			output_norm_select[i] = 0;
	}
	//combine neighbour points
	bool  flag = 0;
	float localMax = 0;
	int   localMaxPos = 0;
	for (int i=0;i<sigLength;i++)
	{
		if(output_norm_select[i]>0)//means candicate encountered
		{
			if (0==flag)
			{
				flag = 1;
				localMax = output_norm_select[i];
				localMaxPos = i;
			}
			else//means flag equals 1
			{
				if (i-localMaxPos>posThres)
				{
					localMax = output_norm_select[i];
					localMaxPos = i;
				}
				else//means combination is needed
				{
					if (output_norm_select[i]<localMax)
					{
						output_norm_select[i]=0;
					}
					else//means localmax need updating
					{
						output_norm_select[localMaxPos] = 0;
						localMax = output_norm_select[i];
						localMaxPos = i;
					}
				}
			}
		}

	}
	//obtain localMaxPos
	int num=0;
	for (int i=0;i<sigLength;i++)
	{
		if(output_norm_select[i]>0)
			num++;
	}
	float *leapPositions = (float*)malloc(num * sizeof(float));
	num = 0;
	for (int i=0;i<sigLength;i++)
	{
		if (output_norm_select[i]>0)
		{
			leapPositions[num++] = i/float(SAMPLING_FREQENCY);
		}
	}
	len =num;
	/*printf("There are %d possible jumping positions laocated at:\n");
	for (int i=0;i<num;i++)
	{
		printf("%.2fs ",leapPositions[i]);
	}
	printf("\n");*/
	


	delete output_norm;
	delete output_norm_select;

    return leapPositions;


	/*FILE* ofp;
	if ((ofp = fopen(argv[2], "w")) == NULL)
	{
	printf("Cannot open output file!\n");
	exit(0);
	}
	for(int n=0; n<sigLength; n++)
	fprintf(ofp, "%f\n", output[n]);
	fclose(ofp);*/
}