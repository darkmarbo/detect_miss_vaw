#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define PI (3.1415926535897932384626433832795)
#define NEIB 5
#define THRESHOLD_44100 0.85
#define THRESHOLD_32000 0.6
#define THRESHOLD_22050 0.6
#define THRESHOLD_16000 0.5
#define THRESHOLD_11025 0.2
#define THRESHOLD_8000  0.2

#define posThres 5 //13.6ms under 44.1kHz

//#define SAMPLING_FREQENCY 44100
//#define MAX_SIG_LENGTH (SAMPLING_FREQENCY*10)

float Input[441000];
float output[441000];
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

int ReadFile(char *wfile,float *freq)
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
			if(head.samplerate>=48000||head.samplerate<0)
			{
				fclose(fp);
				exit(-1);
			}
			SAMFREQ = head.samplerate;
			channel_num = head.channels;
			*freq = SAMFREQ;
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

void fft(float *inputR, float *inputI, int N, float direct)
{
	long sigL, i, j, k, n, period, twoPeriod;
	float tmpR, tmpI, uR, uI, wR, wI;

	sigL = long(pow(float(2), float(N)));

	j = 1;
	for(i=1; i<sigL; i++)
	{
		if(i < j)
		{
			tmpR = inputR[j-1];
			tmpI = inputI[j-1];

			inputR[j-1] = inputR[i-1];
			inputI[j-1] = inputI[i-1];

			inputR[i-1] = tmpR;
			inputI[i-1] = tmpI;
		}

		k = sigL/2;
		while (k < j){ j -=  k;	k /= 2;	}
		j += k;
	}

	for(n=1; n<=N; n++ )
	{  
		twoPeriod = long(pow(float(2), float(n)));
		period = twoPeriod/2;
		uR = 1.0; 
		uI = 0.0; 
		wR = cos( PI/period ); 
		wI = -1.0 * sin( PI/period * direct);

		for(j=0; j<period; j++ ) 
		{  
			for(i=j; i<sigL; i+=twoPeriod)
			{
				tmpR = inputR[i+period]*uR - inputI[i+period]*uI;
				tmpI = inputR[i+period]*uI + inputI[i+period]*uR;

				inputR[i+period] = inputR[i] - tmpR; 
				inputI[i+period] = inputI[i] - tmpI; 
				inputR[i] += tmpR ; 
				inputI[i] += tmpI; 
			}
			tmpR = uR*wR - uI*wI; 
			tmpI = uR*wI + uI*wR; 
			uR = tmpR; 
			uI = tmpI; 
		} 
	} 
}

int combine(float **specGram,int numFrame,int numChan,int numIter)
{
	int comChan = numChan;
	for (int i=0;i<numIter;i++)
	{
		comChan = ceil(float(comChan)/float(2));
		for (int n = 0;n<numFrame;n++)
		{
			if (comChan%2==0)
			{
				for (int j=0;j<comChan;j++)
				{
					specGram[n][j] = specGram[n][2*j]+specGram[n][2*j+1];
				}
			}
			else
			{
				for (int j=0;j<comChan-1;j++)
				{
					specGram[n][j] = specGram[n][2*j]+specGram[n][2*j+1];
				}
				specGram[n][comChan-1] = 2*specGram[n][2*(comChan-1)+1];
			}
		}
	}
	return comChan;
}

float * leapDetect(char * filename,int *len);
void main(int argc, char *argv[])
{
	char* filename = argv[1];
	int len;
	float* leapPositions = leapDetect(filename,&len);
	if (leapPositions==NULL)
	{
		//printf("read error...");
		//return;
		len = 0;
	}

	FILE* fp;
	if((fp=fopen(argv[2],"a+"))==NULL)
	{
		printf("Please provide the second parameter...");
		return;
	}
	else
	{
		fprintf(fp,"\n%d possible jumping positions:\n",len);
		if(len>0)
		{
			for (int i=0;i<len;i++)
			{			
				fprintf(fp,"%.2fs ",leapPositions[2*i]);
			}
			fprintf(fp,"\n");
			fprintf(fp,"%d possible jumping possibilities:\n",len);
			for (int i=0;i<len;i++)
			{			
				fprintf(fp,"%.2f  ",leapPositions[2*i+1]);
			}
			fprintf(fp,"\n");
		}
		
	}

	/*printf("There are %d possible jumping positions laocated at:\n",len);
	for (int i=0;i<len;i++)
	{
	printf("%.2fs ",leapPositions[i]);
	}
	printf("\n");
	*/



}

float * leapDetect(char * filename,int *len)
{
	float *filter;
	float beta;
	int fLength;

	//读入语音数据
	float freq;
	sigLength = ReadFile(filename,&freq);
	if (sigLength<1)
	{
		return NULL;
	}
	for (int i=0;i<sigLength;i++)
	{
		Input[i] = allbuf[i]/float(32768);
	}
	//frame+window
	int tmp = freq*0.02;
	int wlen = tmp - tmp%4;
	float *win = (float*)malloc(wlen * sizeof(float));
	for (int i=0;i<wlen;i++)
	{
		win[i] = sin((0.5+i)/float(wlen)*PI);
	}
	int nfram = floor(float(sigLength)/float(wlen)*float(2));
	int swinLen = ceil(float((nfram+1)*wlen/2));
	float *swin= (float*)malloc(swinLen * sizeof(float));
	for (int k=0;k<swinLen;k++)
	{
		swin[k] = 0;
	}
	//swin(t*wlen/2+1:t*wlen/2+wlen)=swin(t*wlen/2+1:t*wlen/2+wlen)+win.^2;
	for(int t=0;t<nfram;t++)
	{
		for (int k=t*wlen/2;k<t*wlen/2+wlen;k++)
		{
			swin[k] = swin[k]+pow(float(win[k-t*wlen/2]),float(2));
		}
	}
	//短时傅里叶变换(STFT)
	float *buffer = (float*)malloc(wlen * sizeof(float));
	float *inputR, *inputI;
	float **specGram;
	int N = int(ceil(log(float(wlen))/log(2.0)));
	long sigL = long(pow(float(2), float(N)));
	inputR = (float*)malloc(sigL * sizeof(float));
	inputI = (float*)malloc(sigL * sizeof(float));
	int numChan = ceil(pow(float(2),float(N))/float(2))+1;
	specGram = (float**)malloc(nfram * sizeof(float*));
	float** specGram_22p05 = (float**)malloc(nfram * sizeof(float*));
	for (int i=0;i<nfram;i++)
	{
		specGram[i] = (float *)malloc(numChan*sizeof(float));
		specGram_22p05[i] = (float *)malloc(numChan*sizeof(float));
	}
	for (int i=0;i<nfram;i++)
	{
		for (int j=0;j<wlen;j++)
		{
			buffer[j] = 0;
		}
		for (int j=0;j<wlen;j++)
		{
			buffer[j] = win[j]*Input[i*wlen/2+j]/float(wlen)/swin[i*wlen/2+j];
		}
		for(int n=0; n<sigL; n++)
		{
			inputR[n] = (n < wlen) ? buffer[n] : 0;
			inputI[n] = 0;
		}
		fft(inputR, inputI, N, 1);

		for(int n=0; n<numChan; n++)
			specGram[i][n] = sqrt(inputR[n]*inputR[n]+inputI[n]*inputI[n]) + 1e-30;
	}
	//channel combination,通道合并，减少特征维数
	int comChan;
	if (freq>=32000)
	{
		comChan = combine(specGram,nfram,numChan,6);
	}
	else if(freq>=16000)
	{
		comChan = combine(specGram,nfram,numChan,5);
	}
	else
	{
		comChan = combine(specGram,nfram,numChan,4);
	}
	for (int i=0;i<nfram;i++)
	{
		for (int j=0;j<comChan;j++)
		{
			specGram[i][j] = log(specGram[i][j]);
			specGram_22p05[i][j] = specGram[i][j];
		}
	}
	for (int i=NEIB+1;i<nfram-NEIB;i++)//第一帧的fft过大，从第二帧开始
	{
		float tmp,tmp_pre,tmp_pos;
		for (int k=0;k<comChan;k++)
		{
			tmp = 0;tmp_pre = 0;tmp_pos = 0;
			for (int m=-NEIB;m<=NEIB;m++)
			{
				tmp += (specGram[i+m][k]);
			}
			tmp /= float(NEIB*2+1);
			specGram[i][comChan+k] = specGram[i][k] - tmp;
			for (int m=-NEIB;m<=0;m++)
			{
				tmp_pre += (specGram[i+m][k]);
			}
			tmp_pre /= float(NEIB+1);
			for (int m=0;m<=NEIB;m++)
			{
				tmp_pos += (specGram[i+m][k]);
			}
			tmp_pos /= float(NEIB+1);
			float tmp_max = tmp;
			if (tmp_pre>tmp_max) tmp_max = tmp_pre;
			if (tmp_pos>tmp_max) tmp_max = tmp_pos;
			specGram_22p05[i][comChan+k] = specGram[i][k] - tmp_max;
		}
	}
	for(int i=0;i<NEIB+1;i++)
	{
		for (int k=0;k<comChan;k++)
		{
			specGram[i][comChan+k] = specGram[NEIB+1][comChan+k];
			specGram_22p05[i][comChan+k] = specGram_22p05[NEIB+1][comChan+k];
		}
	}
	for(int i=nfram-NEIB;i<nfram;i++)
	{
		for (int k=0;k<comChan;k++)
		{
			specGram[i][comChan+k] = specGram[nfram-NEIB-1][comChan+k];
			specGram_22p05[i][comChan+k] = specGram_22p05[nfram-NEIB-1][comChan+k];
		}
	}
	if (44100==freq)
	{
		//detectByNet
		float* mn = (float*)malloc(comChan*2*sizeof(float));
		mn[0] = -4.05453294545455; mn[1] = -5.48252085454546; mn[2] = -5.72797300606060; mn[3] = -6.16047056363637; mn[4] = -6.74085147878787;
		mn[5] = -7.78524236363636; mn[6] = -8.55422266060606; mn[7] = -9.02388830303030; mn[8] = -7.42303760606061;
		mn[9]= 0.180912448484848; mn[10]= 0.336476709090909; mn[11]= 0.466350854545455; mn[12]= 0.531009539393940; mn[13]= 0.540033933333333;
		mn[14]= 0.647626660606061; mn[15]= 0.670444054545454; mn[16]= 0.709887721212121; mn[17]= 0.532052121212121;

		float* sd = (float*)malloc(comChan*2*sizeof(float));
		sd[0] = 2.31006062520455; sd[1] = 2.47835150556818; sd[2] = 2.71588831528287; sd[3] = 2.57765743144463; sd[4] = 2.32985095519208;
		sd[5] = 2.03524852955554; sd[6] = 1.94477913339129; sd[7] = 1.93682216282548; sd[8] = 2.32976456969904;
		sd[9]= 0.580552064209108; sd[10]= 0.773567176989994; sd[11]= 0.830032413421823; sd[12]= 0.844603341771980; sd[13]= 0.874331791461113;
		sd[14]= 1.07674535273000; sd[15]= 1.24198973438201; sd[16]= 1.35821088416416; sd[17]= 0.930090214744849; 

		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = (specGram[i][j]-mn[j])/sd[j];
			}
		}
		float w1[7][19];
		float w2[8];

		w1[0][0] = 1.659890; w1[0][1] = -1.899934; w1[0][2] = -1.302368; w1[0][3] = -0.590487; w1[0][4] = -0.897917; w1[0][5] = -0.315630; w1[0][6] = 0.718890; w1[0][7] = -1.218169; w1[0][8] = -0.459738; w1[0][9] = 0.567702; w1[0][10] = 0.891833; w1[0][11] = -0.278441; w1[0][12] = -1.319041; w1[0][13] = -1.946008; w1[0][14] = -1.728517; w1[0][15] = -1.074936; w1[0][16] = 1.112241; w1[0][17] = -1.740831; w1[0][18] = 1.896380;
		w1[1][0] = -1.036665; w1[1][1] = 1.120828; w1[1][2] = -0.435549; w1[1][3] = 1.254150; w1[1][4] = -1.057148; w1[1][5] = 1.436789; w1[1][6] = -0.473021; w1[1][7] = -1.371013; w1[1][8] = 1.151134; w1[1][9] = 1.334084; w1[1][10] = -1.249728; w1[1][11] = 0.594796; w1[1][12] = 1.846145; w1[1][13] = 1.609601; w1[1][14] = 1.862827; w1[1][15] = 0.255387; w1[1][16] = 0.530824; w1[1][17] = -0.520996; w1[1][18] = 1.116251;
		w1[2][0] = -0.327501; w1[2][1] = -0.272155; w1[2][2] = 1.652271; w1[2][3] = -1.955790; w1[2][4] = -2.113324; w1[2][5] = -0.089608; w1[2][6] = 1.255358; w1[2][7] = 0.976108; w1[2][8] = 1.964346; w1[2][9] = 0.427717; w1[2][10] = 0.773066; w1[2][11] = 0.209933; w1[2][12] = -0.931418; w1[2][13] = -1.229129; w1[2][14] = 0.322151; w1[2][15] = 1.315760; w1[2][16] = 0.291012; w1[2][17] = -0.207318; w1[2][18] = -0.304712;
		w1[3][0] = 0.383179; w1[3][1] = 0.188846; w1[3][2] = 1.540378; w1[3][3] = -1.181548; w1[3][4] = 1.476735; w1[3][5] = -0.716458; w1[3][6] = 0.074650; w1[3][7] = -1.254154; w1[3][8] = -0.711887; w1[3][9] = -0.134310; w1[3][10] = 1.065821; w1[3][11] = 1.387420; w1[3][12] = -0.967551; w1[3][13] = -0.707154; w1[3][14] = -0.370813; w1[3][15] = -0.527951; w1[3][16] = -0.584072; w1[3][17] = 0.375888; w1[3][18] = -0.046930;
		w1[4][0] = 1.683297; w1[4][1] = 1.132011; w1[4][2] = -0.455457; w1[4][3] = -1.006809; w1[4][4] = 0.359972; w1[4][5] = 1.206223; w1[4][6] = 0.982448; w1[4][7] = -0.434168; w1[4][8] = -1.144199; w1[4][9] = -0.223315; w1[4][10] = -0.640076; w1[4][11] = 0.955431; w1[4][12] = 0.664654; w1[4][13] = 0.175582; w1[4][14] = -0.496320; w1[4][15] = -0.270694; w1[4][16] = 0.125251; w1[4][17] = -3.081091; w1[4][18] = -1.468346;
		w1[5][0] = 0.682421; w1[5][1] = 1.621957; w1[5][2] = 1.257140; w1[5][3] = -1.601064; w1[5][4] = 1.057978; w1[5][5] = 0.810792; w1[5][6] = -0.729705; w1[5][7] = -1.440738; w1[5][8] = 2.054717; w1[5][9] = -1.123714; w1[5][10] = 0.616779; w1[5][11] = 0.218053; w1[5][12] = 1.117283; w1[5][13] = -0.950200; w1[5][14] = -0.533759; w1[5][15] = -0.936232; w1[5][16] = -1.333307; w1[5][17] = 1.718262; w1[5][18] = 0.659722;
		w1[6][0] = -0.584653; w1[6][1] = 0.767554; w1[6][2] = -1.512569; w1[6][3] = 0.757979; w1[6][4] = 0.014799; w1[6][5] = -1.719427; w1[6][6] = -0.570267; w1[6][7] = 1.168406; w1[6][8] = 1.801403; w1[6][9] = 0.274059; w1[6][10] = -0.309939; w1[6][11] = 1.668894; w1[6][12] = -1.104416; w1[6][13] = 0.139654; w1[6][14] = 0.032847; w1[6][15] = -0.196317; w1[6][16] = 0.759955; w1[6][17] = 2.268718; w1[6][18] = -0.492387;


		w2[0] = 0.747940; w2[1] = -0.210797; w2[2] = -1.296694; w2[3] = 1.567248; w2[4] = -2.821456; w2[5] = -2.748451; w2[6] = 1.707678; w2[7] = 2.142546; 


		float *score = (float *)malloc(nfram*sizeof(float));
		for (int n=0;n<nfram;n++)
		{
			float hid[7];
			for (int i=0;i<7;i++) {hid[i]=0;}
			for (int i=0;i<7;i++)
			{
				for (int j=0;j<19;j++)
				{
					if (0==j) hid[i] += w1[i][j];
					else hid[i] += specGram[n][j-1]*w1[i][j];
				}
				hid[i] = 1/(1+exp(-hid[i]));
			}
			float out=0;
			for (int i=0;i<8;i++)
			{
				if(0==i) out += w2[i];
				else out += hid[i-1]*w2[i];
			}
			out = 1/(1+exp(-out));
			score[n] = out;
		}

		float* output_norm_select = (float*)malloc(nfram * sizeof(float));
		//pick out candidate points
		for (int i=0;i<nfram;i++)
		{
			if (score[i]> THRESHOLD_44100)
			{
				output_norm_select[i] = score[i];
			}
			else
				output_norm_select[i] = 0;
		}
		//合并相邻的帧
		int flag = 0;
		float localMax = 0;
		int localPos = 0;
		for (int i=0;i<nfram;i++)
		{
			if (output_norm_select[i]>0)
			{
				if (flag>0)
				{
					if (output_norm_select[i]>localMax)
					{
						output_norm_select[localPos] = 0;
						localPos = i;
						localMax = output_norm_select[i];
					}
					else
					{
						output_norm_select[i] = 0;
					}
				}
				else
				{
					flag = 1;
					localPos = i;
					localMax = output_norm_select[i];
				}
			}
			else
			{
				flag = 0;
			}
		}

		//obtain localMaxPos
		int num=0;
		for (int i=0;i<nfram;i++)
		{
			if(output_norm_select[i]>0)
				num++;
		}
		if (num>0)
		{
			float *leapPositions = (float*)malloc(num*2 * sizeof(float));
			num = 0;
			for (int i=0;i<nfram;i++)
			{
				if (output_norm_select[i]>0)
				{
					leapPositions[2*num] = i*0.01;
					leapPositions[2*num+1] = output_norm_select[i];
					num++;
				}
			}
			*len = num;
			return leapPositions;
		}
		else
		{
			float *leapPositions = NULL;
			*len = num;
			return leapPositions;
		}
	}
	if (32000==freq)
	{
		//detectByNet
		float* mn = (float*)malloc(comChan*2*sizeof(float));
		mn[0] = -3.549642; mn[1] = -4.847797; mn[2] = -5.679484; mn[3] = -5.815553; mn[4] = -6.126271;
		mn[5] = -6.326262; mn[6] = -7.183693; mn[7] = -8.163559; mn[8] = -6.688313;
		mn[9] = 0.115492; mn[10] = 0.212633; mn[11] = 0.295018; mn[12] = 0.343536; mn[13] = 0.387155;
		mn[14] = 0.366545; mn[15] = 0.403909; mn[16] = 0.424791; mn[17] = 0.386625; 

		float* sd = (float*)malloc(comChan*2*sizeof(float));
		sd[0] = 2.280795; sd[1] = 2.426023; sd[2] = 2.371193; sd[3] = 2.519315; sd[4] = 2.396306;
		sd[5] = 2.271404; sd[6] = 2.003222; sd[7] = 1.815204; sd[8] = 2.449061;
		sd[9] = 0.610538; sd[10] = 0.881275; sd[11] = 0.869386; sd[12] = 0.910046; sd[13] = 0.888033;
		sd[14] = 0.887948; sd[15] = 0.871007; sd[16] = 0.987936; sd[17] = 0.896978; 

		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = (specGram[i][j]-mn[j])/sd[j];
			}
		}
		float w1[7][19];
		float w2[8];
		w1[0][0] = 1.554315; w1[0][1] = -1.944684; w1[0][2] = -1.369835; w1[0][3] = -0.641385; w1[0][4] = -0.933820; w1[0][5] = -0.309842; w1[0][6] = 0.714271; w1[0][7] = -1.161506; w1[0][8] = -0.459674; w1[0][9] = 0.573807; w1[0][10] = 0.768666; w1[0][11] = -0.385442; w1[0][12] = -1.396251; w1[0][13] = -1.993578; w1[0][14] = -1.744430; w1[0][15] = -1.103591; w1[0][16] = 1.144320; w1[0][17] = -1.746152; w1[0][18] = 1.881424;
		w1[1][0] = -1.057904; w1[1][1] = 0.978862; w1[1][2] = -0.502669; w1[1][3] = 1.273702; w1[1][4] = -1.044211; w1[1][5] = 1.477084; w1[1][6] = -0.270468; w1[1][7] = -1.130724; w1[1][8] = 1.164752; w1[1][9] = 1.327037; w1[1][10] = -1.108587; w1[1][11] = 0.690755; w1[1][12] = 1.946029; w1[1][13] = 1.586987; w1[1][14] = 1.930221; w1[1][15] = 0.355006; w1[1][16] = 0.538249; w1[1][17] = -0.788343; w1[1][18] = 1.105764;
		w1[2][0] = -0.125759; w1[2][1] = -0.216186; w1[2][2] = 1.553678; w1[2][3] = -2.036580; w1[2][4] = -1.987120; w1[2][5] = 0.000241; w1[2][6] = 1.265174; w1[2][7] = 0.677408; w1[2][8] = 1.554764; w1[2][9] = 0.558171; w1[2][10] = 1.469529; w1[2][11] = 0.573577; w1[2][12] = -0.700225; w1[2][13] = -0.939980; w1[2][14] = 0.482161; w1[2][15] = 1.249000; w1[2][16] = 0.138473; w1[2][17] = -0.373492; w1[2][18] = -0.061005;
		w1[3][0] = 0.344890; w1[3][1] = -0.035144; w1[3][2] = 1.381583; w1[3][3] = -1.655507; w1[3][4] = 0.669029; w1[3][5] = -1.605578; w1[3][6] = -0.341659; w1[3][7] = -1.549338; w1[3][8] = -1.235892; w1[3][9] = -1.138976; w1[3][10] = -0.630264; w1[3][11] = 0.573900; w1[3][12] = -1.211799; w1[3][13] = -1.096842; w1[3][14] = -0.946913; w1[3][15] = -0.481404; w1[3][16] = -0.912137; w1[3][17] = -0.253788; w1[3][18] = -0.772514;
		w1[4][0] = 1.595257; w1[4][1] = 0.749519; w1[4][2] = -0.244771; w1[4][3] = -1.068207; w1[4][4] = -0.061901; w1[4][5] = 1.074719; w1[4][6] = 1.054635; w1[4][7] = 0.358927; w1[4][8] = -0.857869; w1[4][9] = -0.493528; w1[4][10] = -0.667258; w1[4][11] = 1.531737; w1[4][12] = 0.970784; w1[4][13] = 0.144293; w1[4][14] = 0.044000; w1[4][15] = 0.339916; w1[4][16] = 0.548654; w1[4][17] = -3.590119; w1[4][18] = -1.195748;
		w1[5][0] = 0.799749; w1[5][1] = 1.734394; w1[5][2] = 1.199304; w1[5][3] = -1.458825; w1[5][4] = 1.434467; w1[5][5] = 0.854830; w1[5][6] = -0.556516; w1[5][7] = -1.801654; w1[5][8] = 1.668833; w1[5][9] = -1.079229; w1[5][10] = 1.053146; w1[5][11] = 0.208674; w1[5][12] = 1.488388; w1[5][13] = -0.199707; w1[5][14] = -0.409991; w1[5][15] = -0.755570; w1[5][16] = -1.304539; w1[5][17] = 1.876667; w1[5][18] = 0.759829;
		w1[6][0] = -0.704547; w1[6][1] = 1.069538; w1[6][2] = -1.711291; w1[6][3] = 0.466487; w1[6][4] = -0.071107; w1[6][5] = -1.943429; w1[6][6] = -0.859962; w1[6][7] = 0.637342; w1[6][8] = 1.390368; w1[6][9] = 0.125421; w1[6][10] = 0.034758; w1[6][11] = 1.215409; w1[6][12] = -1.531021; w1[6][13] = -0.156914; w1[6][14] = -0.443926; w1[6][15] = -0.683010; w1[6][16] = 0.189822; w1[6][17] = 2.047275; w1[6][18] = -0.779126;


		w2[0] = 0.318144; w2[1] = -0.614089; w2[2] = -1.368598; w2[3] = 0.052299; w2[4] = -3.262472; w2[5] = -2.849227; w2[6] = 1.360852; w2[7] = 1.605819; 

		float *score = (float *)malloc(nfram*sizeof(float));
		for (int n=0;n<nfram;n++)
		{
			float hid[7];
			for (int i=0;i<7;i++) {hid[i]=0;}
			for (int i=0;i<7;i++)
			{
				for (int j=0;j<19;j++)
				{
					if (0==j) hid[i] += w1[i][j];
					else hid[i] += specGram[n][j-1]*w1[i][j];
				}
				hid[i] = 1/(1+exp(-hid[i]));
			}
			float out=0;
			for (int i=0;i<8;i++)
			{
				if(0==i) out += w2[i];
				else out += hid[i-1]*w2[i];
			}
			out = 1/(1+exp(-out));
			score[n] = out;
		}
		//select candidata points
		float* output_norm_select = (float*)malloc(nfram * sizeof(float));
		//pick out candidate points
		for (int i=0;i<nfram;i++)
		{
			if (score[i]> THRESHOLD_32000)
			{
				output_norm_select[i] = score[i];
			}
			else
				output_norm_select[i] = 0;
		}
		//合并相邻的帧
		int flag = 0;
		float localMax = 0;
		int localPos = 0;
		for (int i=0;i<nfram;i++)
		{
			if (output_norm_select[i]>0)
			{
				if (flag>0)
				{
					if (output_norm_select[i]>localMax)
					{
						output_norm_select[localPos] = 0;
						localPos = i;
						localMax = output_norm_select[i];
					}
					else
					{
						output_norm_select[i] = 0;
					}
				}
				else
				{
					flag = 1;
					localPos = i;
					localMax = output_norm_select[i];
				}
			}
			else
			{
				flag = 0;
			}
		}
		//obtain localMaxPos
		int num=0;
		for (int i=0;i<nfram;i++)
		{
			if(output_norm_select[i]>0)
				num++;
		}
		if (num>0)
		{
			float *leapPositions = (float*)malloc(2*num * sizeof(float));
			num = 0;
			for (int i=0;i<nfram;i++)
			{
				if (output_norm_select[i]>0)
				{
					leapPositions[2*num] = i*0.01;
					leapPositions[2*num+1] = output_norm_select[i];
					num++;
				}
				
			}
			*len = num;
			return leapPositions;
		}
		else
		{
			*len = num;
			float *leapPositions = NULL;
			return leapPositions;
		}
	}
	if (22050==freq)
	{
		//detectByNet
		float* mn = (float*)malloc(comChan*2*sizeof(float));
		mn[0] = -3.902824; mn[1] = -4.848543; mn[2] = -5.556586; mn[3] = -6.480794; mn[4] = -6.450456;
		mn[5] = -6.679172; mn[6] = -6.908638; mn[7] = -7.306495; mn[8] = -7.064706;
		mn[9] = -0.365937; mn[10] = -0.374813; mn[11] = -0.320477; mn[12] = -0.205396; mn[13] = -0.226712;
		mn[14] = -0.214537; mn[15] = -0.164160; mn[16] = -0.156791; mn[17] = -0.259914; 

		float* sd = (float*)malloc(comChan*2*sizeof(float));
		sd[0] = 2.395000; sd[1] = 2.355090; sd[2] = 2.499537; sd[3] = 2.230427; sd[4] = 2.407242;
		sd[5] = 2.353789; sd[6] = 2.235740; sd[7] = 2.131429; sd[8] = 2.426038;
		sd[9] = 0.682916; sd[10] = 0.794873; sd[11] = 0.788831; sd[12] = 0.802726; sd[13] = 0.824534;
		sd[14] = 0.872025; sd[15] = 0.871683; sd[16] = 0.838339; sd[17] = 0.889591; 

		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = specGram_22p05[i][j];
			}
		}
		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = (specGram[i][j]-mn[j])/sd[j];
			}
		}
		
		float w1[7][19]; float w2[8];
		w1[0][0] = 1.850193; w1[0][1] = -1.800741; w1[0][2] = -1.168207; w1[0][3] = -0.470150; w1[0][4] = -0.859134; w1[0][5] = -0.244734; w1[0][6] = 0.773409; w1[0][7] = -1.177823; w1[0][8] = -0.421173; w1[0][9] = 0.670767; w1[0][10] = 0.934078; w1[0][11] = -0.258311; w1[0][12] = -1.277129; w1[0][13] = -1.934593; w1[0][14] = -1.717255; w1[0][15] = -1.036222; w1[0][16] = 1.135520; w1[0][17] = -1.727993; w1[0][18] = 1.922468;
		w1[1][0] = -0.829232; w1[1][1] = 0.806976; w1[1][2] = -0.564202; w1[1][3] = 1.238317; w1[1][4] = -1.308358; w1[1][5] = 1.311185; w1[1][6] = -0.300646; w1[1][7] = -1.086142; w1[1][8] = 1.517808; w1[1][9] = 1.268332; w1[1][10] = -0.855009; w1[1][11] = 0.848539; w1[1][12] = 2.054103; w1[1][13] = 1.768820; w1[1][14] = 1.924363; w1[1][15] = 0.378332; w1[1][16] = 0.686123; w1[1][17] = -0.358010; w1[1][18] = 1.264623;
		w1[2][0] = -0.077845; w1[2][1] = -0.009015; w1[2][2] = 1.802392; w1[2][3] = -1.894917; w1[2][4] = -1.697259; w1[2][5] = 0.243222; w1[2][6] = 1.408277; w1[2][7] = 0.843166; w1[2][8] = 1.460291; w1[2][9] = 0.721871; w1[2][10] = 1.494942; w1[2][11] = 0.697915; w1[2][12] = -0.719033; w1[2][13] = -0.811084; w1[2][14] = 0.619839; w1[2][15] = 1.246393; w1[2][16] = 0.030052; w1[2][17] = -0.609619; w1[2][18] = 0.033147;
		w1[3][0] = 0.990657; w1[3][1] = -0.058095; w1[3][2] = 2.058387; w1[3][3] = -0.848902; w1[3][4] = 0.280068; w1[3][5] = -1.506338; w1[3][6] = -0.400359; w1[3][7] = -1.332186; w1[3][8] = 0.160593; w1[3][9] = -0.616183; w1[3][10] = 0.152810; w1[3][11] = 0.906941; w1[3][12] = -0.673347; w1[3][13] = -0.928999; w1[3][14] = -0.574187; w1[3][15] = -0.820168; w1[3][16] = -1.033410; w1[3][17] = 0.354397; w1[3][18] = 0.147659;
		w1[4][0] = 2.516375; w1[4][1] = -0.100682; w1[4][2] = -0.445132; w1[4][3] = -0.821085; w1[4][4] = -0.624005; w1[4][5] = 0.690426; w1[4][6] = 0.754878; w1[4][7] = -0.179494; w1[4][8] = -0.232273; w1[4][9] = -0.477097; w1[4][10] = -0.522696; w1[4][11] = 2.407863; w1[4][12] = 2.095945; w1[4][13] = 0.477032; w1[4][14] = 0.183144; w1[4][15] = -0.243676; w1[4][16] = 0.152195; w1[4][17] = -2.632502; w1[4][18] = -0.459655;
		w1[5][0] = 0.017347; w1[5][1] = 2.021557; w1[5][2] = 1.600916; w1[5][3] = -1.539583; w1[5][4] = 1.247627; w1[5][5] = 0.802738; w1[5][6] = -1.058523; w1[5][7] = -1.913920; w1[5][8] = 1.182551; w1[5][9] = -1.235950; w1[5][10] = 1.256677; w1[5][11] = 0.212988; w1[5][12] = 1.100045; w1[5][13] = -0.249218; w1[5][14] = -0.073094; w1[5][15] = -0.762533; w1[5][16] = -1.221449; w1[5][17] = 1.500106; w1[5][18] = 1.002745;
		w1[6][0] = -0.725541; w1[6][1] = 1.357306; w1[6][2] = -1.460795; w1[6][3] = 0.474694; w1[6][4] = 0.104207; w1[6][5] = -1.770356; w1[6][6] = -0.681441; w1[6][7] = 0.921660; w1[6][8] = 1.369120; w1[6][9] = 0.171520; w1[6][10] = 0.080981; w1[6][11] = 1.468418; w1[6][12] = -1.513892; w1[6][13] = -0.044729; w1[6][14] = -0.251888; w1[6][15] = -0.406587; w1[6][16] = 0.522655; w1[6][17] = 1.824796; w1[6][18] = -0.757678;


		w2[0] = -0.017560; w2[1] = -0.906432; w2[2] = -1.864789; w2[3] = 0.578560; w2[4] = -3.344397; w2[5] = -3.152771; w2[6] = 2.510419; w2[7] = 0.913048; 

		float *score = (float *)malloc(nfram*sizeof(float));
		for (int n=0;n<nfram;n++)
		{
			float hid[7];
			for (int i=0;i<7;i++) {hid[i]=0;}
			for (int i=0;i<7;i++)
			{
				for (int j=0;j<19;j++)
				{
					if (0==j) hid[i] += w1[i][j];
					else hid[i] += specGram[n][j-1]*w1[i][j];
				}
				hid[i] = 1/(1+exp(-hid[i]));
			}
			float out=0;
			for (int i=0;i<8;i++)
			{
				if(0==i) out += w2[i];
				else out += hid[i-1]*w2[i];
			}
			out = 1/(1+exp(-out));
			score[n] = out;
		}
		/*FILE *fpout = fopen("prob.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpout,"%f\n",score[n]);
		}
		fclose(fpout);*/
		//select candidata points
		float* output_norm_select = (float*)malloc(nfram * sizeof(float));
		//pick out candidate points
		for (int i=0;i<nfram;i++)
		{
			if (score[i]> THRESHOLD_22050)
			{
				output_norm_select[i] = score[i];
			}
			else
				output_norm_select[i] = 0;
		}

		//合并相邻的帧
		int flag = 0;
		float localMax = 0;
		int localPos = 0;
		for (int i=0;i<nfram;i++)
		{
			if (output_norm_select[i]>0)
			{
				if (flag>0)
				{
					if (output_norm_select[i]>localMax)
					{
						output_norm_select[localPos] = 0;
						localPos = i;
						localMax = output_norm_select[i];
					}
					else
					{
						output_norm_select[i] = 0;
					}
				}
				else
				{
					flag = 1;
					localPos = i;
					localMax = output_norm_select[i];
				}
			}
			else
			{
				flag = 0;
			}
		}

		//obtain localMaxPos
		int num=0;
		for (int i=0;i<nfram;i++)
		{
			if(output_norm_select[i]>0)
				num++;
		}
		if (num>0)
		{
			float *leapPositions = (float*)malloc(2*num * sizeof(float));
			num = 0;
			for (int i=0;i<nfram;i++)
			{
				if (output_norm_select[i]>0)
				{
					leapPositions[2*num] = i*0.01;
					leapPositions[2*num+1] = output_norm_select[i];
					num++;
				}
			}
			*len = num;
			return leapPositions;
		}
		else
		{
			*len = num;
			float *leapPositions = NULL;
			return leapPositions;
		}
	}
	if (16000==freq)
	{
		//detectByNet
		float* mn = (float*)malloc(comChan*2*sizeof(float));
		mn[0] = -3.755666; mn[1] = -4.824424; mn[2] = -5.079571; mn[3] = -5.754408; mn[4] = -6.307023;
		mn[5] = -6.560711; mn[6] = -6.513738; mn[7] = -6.863776; mn[8] = -7.321334;
		mn[9] = -0.387941; mn[10] = -0.400048; mn[11] = -0.373688; mn[12] = -0.331960; mn[13] = -0.240749;
		mn[14] = -0.255527; mn[15] = -0.246672; mn[16] = -0.233224; mn[17] = -0.229599; 

		float* sd = (float*)malloc(comChan*2*sizeof(float));
		sd[0] = 2.427287; sd[1] = 2.264693; sd[2] = 2.536748; sd[3] = 2.443007; sd[4] = 2.243614;
		sd[5] = 2.354053; sd[6] = 2.384589; sd[7] = 2.302715; sd[8] = 2.164407;
		sd[9] = 0.702708; sd[10] = 0.794549; sd[11] = 0.883299; sd[12] = 0.785236; sd[13] = 0.811866;
		sd[14] = 0.825035; sd[15] = 0.844339; sd[16] = 0.846162; sd[17] = 0.848053; 

		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = specGram_22p05[i][j];
			}
		}
		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = (specGram[i][j]-mn[j])/sd[j];
			}
		}

		float w1[7][19];
		float w2[8];
		w1[0][0] = 1.717853; w1[0][1] = -1.881388; w1[0][2] = -1.273729; w1[0][3] = -0.572259; w1[0][4] = -0.886289; w1[0][5] = -0.313900; w1[0][6] = 0.697450; w1[0][7] = -1.231288; w1[0][8] = -0.495164; w1[0][9] = 0.577083; w1[0][10] = 0.899408; w1[0][11] = -0.252563; w1[0][12] = -1.293491; w1[0][13] = -1.925980; w1[0][14] = -1.697258; w1[0][15] = -1.085271; w1[0][16] = 1.089601; w1[0][17] = -1.779078; w1[0][18] = 1.925947;
		w1[1][0] = -0.851225; w1[1][1] = 1.275176; w1[1][2] = -0.409696; w1[1][3] = 1.286431; w1[1][4] = -1.112432; w1[1][5] = 1.424424; w1[1][6] = -0.415223; w1[1][7] = -1.111151; w1[1][8] = 1.504401; w1[1][9] = 1.294828; w1[1][10] = -0.761597; w1[1][11] = 0.657971; w1[1][12] = 1.908235; w1[1][13] = 1.537506; w1[1][14] = 1.937876; w1[1][15] = 0.301348; w1[1][16] = 0.699217; w1[1][17] = -0.339485; w1[1][18] = 1.128248;
		w1[2][0] = -0.032719; w1[2][1] = -0.185122; w1[2][2] = 1.758072; w1[2][3] = -1.912575; w1[2][4] = -1.859758; w1[2][5] = 0.176614; w1[2][6] = 1.386388; w1[2][7] = 0.826481; w1[2][8] = 1.510585; w1[2][9] = 0.711069; w1[2][10] = 1.352277; w1[2][11] = 0.566195; w1[2][12] = -0.859109; w1[2][13] = -1.117101; w1[2][14] = 0.423338; w1[2][15] = 1.053157; w1[2][16] = -0.102210; w1[2][17] = -0.652433; w1[2][18] = -0.158101;
		w1[3][0] = 0.749651; w1[3][1] = 0.146973; w1[3][2] = 1.111346; w1[3][3] = -1.445799; w1[3][4] = 1.008068; w1[3][5] = -1.169011; w1[3][6] = -0.148944; w1[3][7] = -0.875640; w1[3][8] = 0.416644; w1[3][9] = -0.684345; w1[3][10] = 0.113276; w1[3][11] = 0.852956; w1[3][12] = -1.132809; w1[3][13] = -1.106936; w1[3][14] = -0.751351; w1[3][15] = -0.279310; w1[3][16] = -0.360456; w1[3][17] = 0.762331; w1[3][18] = -0.523890;
		w1[4][0] = 2.149321; w1[4][1] = 0.997636; w1[4][2] = -0.801129; w1[4][3] = -1.089706; w1[4][4] = -0.013992; w1[4][5] = 0.842331; w1[4][6] = 0.747083; w1[4][7] = 0.168784; w1[4][8] = 0.130138; w1[4][9] = -0.658136; w1[4][10] = -0.996162; w1[4][11] = 1.110897; w1[4][12] = 0.912372; w1[4][13] = 0.262328; w1[4][14] = -0.346564; w1[4][15] = -0.024243; w1[4][16] = 0.471605; w1[4][17] = -2.359111; w1[4][18] = -1.575533;
		w1[5][0] = 0.405831; w1[5][1] = 1.812268; w1[5][2] = 1.640558; w1[5][3] = -1.100304; w1[5][4] = 1.699704; w1[5][5] = 1.174939; w1[5][6] = -0.505529; w1[5][7] = -1.584704; w1[5][8] = 1.309732; w1[5][9] = -0.816573; w1[5][10] = 1.059720; w1[5][11] = 0.247048; w1[5][12] = 1.143340; w1[5][13] = -0.638468; w1[5][14] = -0.476176; w1[5][15] = -1.060633; w1[5][16] = -1.614318; w1[5][17] = 1.299134; w1[5][18] = 0.662491;
		w1[6][0] = -0.310775; w1[6][1] = 0.908844; w1[6][2] = -1.805631; w1[6][3] = 0.404477; w1[6][4] = -0.107741; w1[6][5] = -1.874649; w1[6][6] = -0.796803; w1[6][7] = 0.751181; w1[6][8] = 1.261750; w1[6][9] = 0.134127; w1[6][10] = 0.130478; w1[6][11] = 1.416420; w1[6][12] = -1.485569; w1[6][13] = -0.197933; w1[6][14] = -0.350214; w1[6][15] = -0.538787; w1[6][16] = 0.289983; w1[6][17] = 1.722832; w1[6][18] = -0.801857;


		w2[0] = 0.102381; w2[1] = -0.636818; w2[2] = -1.702851; w2[3] = 0.783791; w2[4] = -2.930482; w2[5] = -2.672489; w2[6] = 1.875465; w2[7] = -0.484247; 

		float *score = (float *)malloc(nfram*sizeof(float));
		for (int n=0;n<nfram;n++)
		{
			float hid[7];
			for (int i=0;i<7;i++) {hid[i]=0;}
			for (int i=0;i<7;i++)
			{
				for (int j=0;j<19;j++)
				{
					if (0==j) hid[i] += w1[i][j];
					else hid[i] += specGram[n][j-1]*w1[i][j];
				}
				hid[i] = 1/(1+exp(-hid[i]));
			}
			float out=0;
			for (int i=0;i<8;i++)
			{
				if(0==i) out += w2[i];
				else out += hid[i-1]*w2[i];
			}
			out = 1/(1+exp(-out));
			score[n] = out;
		}
		/*FILE *fpout = fopen("prob.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpout,"%f\n",score[n]);
		}
		fclose(fpout);*/
		//select candidata points
		float* output_norm_select = (float*)malloc(nfram * sizeof(float));
		//pick out candidate points
		for (int i=0;i<nfram;i++)
		{
			if (score[i]> THRESHOLD_16000)
			{
				output_norm_select[i] = score[i];
			}
			else
				output_norm_select[i] = 0;
		}

		/*FILE *fpSele = fopen("candiSele.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpSele,"%f\n",output_norm_select[n]);
		}
		fclose(fpSele);*/
		//合并相邻的帧
		int flag = 0;
		float localMax = 0;
		int localPos = 0;
		for (int i=0;i<nfram;i++)
		{
			if (output_norm_select[i]>0)
			{
				if (flag>0)
				{
					if (output_norm_select[i]>localMax)
					{
						output_norm_select[localPos] = 0;
						localPos = i;
						localMax = output_norm_select[i];
					}
					else
					{
						output_norm_select[i] = 0;
					}
				}
				else
				{
					flag = 1;
					localPos = i;
					localMax = output_norm_select[i];
				}
			}
			else
			{
				flag = 0;
			}
		}
		//obtain localMaxPos
		int num=0;
		for (int i=0;i<nfram;i++)
		{
			if(output_norm_select[i]>0)
				num++;
		}
		if (num>0)
		{
			float *leapPositions = (float*)malloc(2*num * sizeof(float));
			num = 0;
			for (int i=0;i<nfram;i++)
			{
				if (output_norm_select[i]>0)
				{
					leapPositions[2*num] = i*0.01;
					leapPositions[2*num+1] = output_norm_select[i];
					num++;
				}
			}
			*len = num;
			return leapPositions;
		}
		else
		{
			*len = num;
			float *leapPositions = NULL;
			return leapPositions;
		}
	}
	if (11025==freq)
	{
		//detectByNet
		float* mn = (float*)malloc(comChan*2*sizeof(float));
		mn[0] = -4.296121; mn[1] = -5.372293; mn[2] = -5.543746; mn[3] = -5.809711; mn[4] = -6.134281;
		mn[5] = -6.512179; mn[6] = -7.162570; mn[7] = -7.458425; mn[8] = -6.725479;
		mn[9] = -0.371848; mn[10] = -0.386629; mn[11] = -0.405380; mn[12] = -0.353643; mn[13] = -0.330324;
		mn[14] = -0.304000; mn[15] = -0.206584; mn[16] = -0.216728; mn[17] = -0.340099; 

		float* sd = (float*)malloc(comChan*2*sizeof(float));
		sd[0] = 2.402445; sd[1] = 2.343163; sd[2] = 2.228363; sd[3] = 2.537669; sd[4] = 2.537380;
		sd[5] = 2.408680; sd[6] = 2.181354; sd[7] = 2.266818; sd[8] = 2.579096;
		sd[9] = 0.685498; sd[10] = 0.737665; sd[11] = 0.828645; sd[12] = 0.893027; sd[13] = 0.845139;
		sd[14] = 0.746975; sd[15] = 0.811487; sd[16] = 0.841029; sd[17] = 0.914232; 

		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = specGram_22p05[i][j];
			}
		}
		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = (specGram[i][j]-mn[j])/sd[j];
			}
		}
		/*FILE *fpnd=fopen("normData.txt","w");
		for (int i=0;i<nfram;i++)
		{
		for (int j=0;j<comChan*2;j++)
		{
		fprintf(fpnd,"%f ",specGram[i][j]);
		}
		fprintf(fpnd,"\n");
		}
		fclose(fpnd);*/

		float w1[7][19];
		float w2[8];
		w1[0][0] = 1.699298; w1[0][1] = -1.854584; w1[0][2] = -1.249674; w1[0][3] = -0.568166; w1[0][4] = -0.881972; w1[0][5] = -0.314027; w1[0][6] = 0.694502; w1[0][7] = -1.256982; w1[0][8] = -0.525079; w1[0][9] = 0.575070; w1[0][10] = 0.943447; w1[0][11] = -0.219820; w1[0][12] = -1.289484; w1[0][13] = -1.933201; w1[0][14] = -1.721683; w1[0][15] = -1.109300; w1[0][16] = 1.083211; w1[0][17] = -1.772016; w1[0][18] = 1.919584;
		w1[1][0] = -0.793492; w1[1][1] = 1.232658; w1[1][2] = -0.439390; w1[1][3] = 1.230358; w1[1][4] = -1.075622; w1[1][5] = 1.422116; w1[1][6] = -0.381192; w1[1][7] = -1.057029; w1[1][8] = 1.517128; w1[1][9] = 1.309568; w1[1][10] = -0.874753; w1[1][11] = 0.550773; w1[1][12] = 1.883583; w1[1][13] = 1.654726; w1[1][14] = 1.972487; w1[1][15] = 0.315417; w1[1][16] = 0.841215; w1[1][17] = -0.184919; w1[1][18] = 1.207521;
		w1[2][0] = -0.129766; w1[2][1] = -0.220828; w1[2][2] = 1.725351; w1[2][3] = -1.917952; w1[2][4] = -1.987240; w1[2][5] = 0.017949; w1[2][6] = 1.254716; w1[2][7] = 0.679337; w1[2][8] = 1.362238; w1[2][9] = 0.538859; w1[2][10] = 1.439160; w1[2][11] = 0.599040; w1[2][12] = -0.750309; w1[2][13] = -1.203779; w1[2][14] = 0.211192; w1[2][15] = 1.059821; w1[2][16] = -0.147341; w1[2][17] = -0.779331; w1[2][18] = -0.362756;
		w1[3][0] = 1.375335; w1[3][1] = 0.272708; w1[3][2] = 1.138299; w1[3][3] = -1.867944; w1[3][4] = 1.290829; w1[3][5] = -0.826402; w1[3][6] = 0.106629; w1[3][7] = -0.925156; w1[3][8] = 0.075984; w1[3][9] = -0.150142; w1[3][10] = -0.505924; w1[3][11] = 0.721891; w1[3][12] = -1.407011; w1[3][13] = -0.946669; w1[3][14] = -0.521961; w1[3][15] = -0.692444; w1[3][16] = -0.441295; w1[3][17] = 1.202381; w1[3][18] = 0.013652;
		w1[4][0] = 2.163872; w1[4][1] = 0.875587; w1[4][2] = -0.634467; w1[4][3] = -1.622353; w1[4][4] = -0.232017; w1[4][5] = 0.847636; w1[4][6] = 0.810930; w1[4][7] = -0.012016; w1[4][8] = -0.254308; w1[4][9] = -0.562339; w1[4][10] = -0.766549; w1[4][11] = 1.709835; w1[4][12] = 1.092939; w1[4][13] = 0.363615; w1[4][14] = -0.076082; w1[4][15] = -0.219941; w1[4][16] = 0.609324; w1[4][17] = -2.339063; w1[4][18] = -1.254209;
		w1[5][0] = 0.533055; w1[5][1] = 1.897295; w1[5][2] = 1.697509; w1[5][3] = -1.008232; w1[5][4] = 1.793506; w1[5][5] = 1.327296; w1[5][6] = -0.348238; w1[5][7] = -1.453155; w1[5][8] = 1.521283; w1[5][9] = -0.660505; w1[5][10] = 1.082906; w1[5][11] = 0.255317; w1[5][12] = 1.147420; w1[5][13] = -0.683169; w1[5][14] = -0.709395; w1[5][15] = -1.107992; w1[5][16] = -1.499903; w1[5][17] = 1.255242; w1[5][18] = 0.390259;
		w1[6][0] = -0.454944; w1[6][1] = 1.100056; w1[6][2] = -1.599177; w1[6][3] = 0.566557; w1[6][4] = -0.024781; w1[6][5] = -1.787523; w1[6][6] = -0.649797; w1[6][7] = 0.911533; w1[6][8] = 1.391232; w1[6][9] = 0.218016; w1[6][10] = 0.192963; w1[6][11] = 1.505733; w1[6][12] = -1.348256; w1[6][13] = -0.082173; w1[6][14] = -0.273986; w1[6][15] = -0.404973; w1[6][16] = 0.513871; w1[6][17] = 1.939056; w1[6][18] = -0.711757;


		w2[0] = -0.009254; w2[1] = -0.202160; w2[2] = -2.140949; w2[3] = 0.476329; w2[4] = -3.246720; w2[5] = -2.749925; w2[6] = 1.765833; w2[7] = 0.769840; 

		float *score = (float *)malloc(nfram*sizeof(float));
		for (int n=0;n<nfram;n++)
		{
			float hid[7];
			for (int i=0;i<7;i++) {hid[i]=0;}
			for (int i=0;i<7;i++)
			{
				for (int j=0;j<19;j++)
				{
					if (0==j) hid[i] += w1[i][j];
					else hid[i] += specGram[n][j-1]*w1[i][j];
				}
				hid[i] = 1/(1+exp(-hid[i]));
			}
			float out=0;
			for (int i=0;i<8;i++)
			{
				if(0==i) out += w2[i];
				else out += hid[i-1]*w2[i];
			}
			out = 1/(1+exp(-out));
			score[n] = out;
		}
		/*FILE *fpout = fopen("prob.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpout,"%f\n",score[n]);
		}
		fclose(fpout);*/
		//select candidata points
		float* output_norm_select = (float*)malloc(nfram * sizeof(float));
		//pick out candidate points
		for (int i=0;i<nfram;i++)
		{
			if (score[i]> THRESHOLD_11025)
			{
				output_norm_select[i] = score[i];
			}
			else
				output_norm_select[i] = 0;
		}

		/*FILE *fpSele = fopen("candiSele.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpSele,"%f\n",output_norm_select[n]);
		}
		fclose(fpSele);*/

		//合并相邻的帧
		int flag = 0;
		float localMax = 0;
		int localPos = 0;
		for (int i=0;i<nfram;i++)
		{
			if (output_norm_select[i]>0)
			{
				if (flag>0)
				{
					if (output_norm_select[i]>localMax)
					{
						output_norm_select[localPos] = 0;
						localPos = i;
						localMax = output_norm_select[i];
					}
					else
					{
						output_norm_select[i] = 0;
					}
				}
				else
				{
					flag = 1;
					localPos = i;
					localMax = output_norm_select[i];
				}
			}
			else
			{
				flag = 0;
			}
		}
		//obtain localMaxPos
		int num=0;
		for (int i=0;i<nfram;i++)
		{
			if(output_norm_select[i]>0)
				num++;
		}
		if (num>0)
		{
			float *leapPositions = (float*)malloc(2*num * sizeof(float));
			num = 0;
			for (int i=0;i<nfram;i++)
			{
				if (output_norm_select[i]>0)
				{
					leapPositions[2*num] = i*0.01;
					leapPositions[2*num+1] = output_norm_select[i];
					num++;
				}
			}
			*len = num;
			return leapPositions;
		}
		else
		{
			*len = num;
			float *leapPositions = NULL;
			return leapPositions;
		}
	}
	if (8000==freq)
	{
		//detectByNet
		float* mn = (float*)malloc(comChan*2*sizeof(float));
		mn[0] = -4.294712; mn[1] = -4.848049; mn[2] = -5.636290; mn[3] = -5.549448; mn[4] = -5.789009;
		mn[5] = -5.879256; mn[6] = -6.416440; mn[7] = -6.762676; mn[8] = -6.575821;
		mn[9] = -0.351840; mn[10] = -0.483815; mn[11] = -0.372991; mn[12] = -0.422260; mn[13] = -0.387761;
		mn[14] = -0.359765; mn[15] = -0.348417; mn[16] = -0.324813; mn[17] = -0.371403; 

		float* sd = (float*)malloc(comChan*2*sizeof(float));
		sd[0] = 2.292559; sd[1] = 2.605435; sd[2] = 2.280849; sd[3] = 2.233843; sd[4] = 2.485573;
		sd[5] = 2.581526; sd[6] = 2.450681; sd[7] = 2.381851; sd[8] = 2.383784;
		sd[9] = 0.678562; sd[10] = 0.789331; sd[11] = 0.773599; sd[12] = 0.856154; sd[13] = 0.925415;
		sd[14] = 0.883449; sd[15] = 0.846466; sd[16] = 0.756297; sd[17] = 0.893692; 

		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = specGram_22p05[i][j];
			}
		}
		for (int i=0;i<nfram;i++)
		{
			for (int j=0;j<comChan*2;j++)
			{
				specGram[i][j] = (specGram[i][j]-mn[j])/sd[j];
			}
		}

		float w1[7][19];
		float w2[8];
		w1[0][0] = 1.696977; w1[0][1] = -1.910582; w1[0][2] = -1.289477; w1[0][3] = -0.586082; w1[0][4] = -0.894511; w1[0][5] = -0.327698; w1[0][6] = 0.695519; w1[0][7] = -1.230421; w1[0][8] = -0.492124; w1[0][9] = 0.558542; w1[0][10] = 0.831687; w1[0][11] = -0.231269; w1[0][12] = -1.327782; w1[0][13] = -1.945124; w1[0][14] = -1.725497; w1[0][15] = -1.108722; w1[0][16] = 1.068828; w1[0][17] = -1.770522; w1[0][18] = 1.906657;
		w1[1][0] = -1.077386; w1[1][1] = 1.071774; w1[1][2] = -0.559033; w1[1][3] = 1.023306; w1[1][4] = -1.358982; w1[1][5] = 1.237654; w1[1][6] = -0.581700; w1[1][7] = -1.289204; w1[1][8] = 1.327232; w1[1][9] = 1.134181; w1[1][10] = -0.858611; w1[1][11] = 0.781825; w1[1][12] = 1.861507; w1[1][13] = 1.635529; w1[1][14] = 1.979089; w1[1][15] = 0.421517; w1[1][16] = 0.648531; w1[1][17] = -0.338683; w1[1][18] = 1.243623;
		w1[2][0] = -0.000786; w1[2][1] = 0.000668; w1[2][2] = 1.919093; w1[2][3] = -1.767126; w1[2][4] = -1.652808; w1[2][5] = 0.322304; w1[2][6] = 1.489395; w1[2][7] = 0.864488; w1[2][8] = 1.480485; w1[2][9] = 0.833245; w1[2][10] = 1.431795; w1[2][11] = 0.649860; w1[2][12] = -0.843309; w1[2][13] = -1.078255; w1[2][14] = 0.322039; w1[2][15] = 1.034350; w1[2][16] = -0.091883; w1[2][17] = -0.666897; w1[2][18] = -0.233557;
		w1[3][0] = 1.001712; w1[3][1] = 0.376254; w1[3][2] = 1.616970; w1[3][3] = -1.580954; w1[3][4] = 0.683782; w1[3][5] = -1.210836; w1[3][6] = 0.104348; w1[3][7] = -0.697667; w1[3][8] = 0.667962; w1[3][9] = -0.554258; w1[3][10] = -0.255603; w1[3][11] = 0.673624; w1[3][12] = -1.326966; w1[3][13] = -1.092857; w1[3][14] = -0.580821; w1[3][15] = -0.308932; w1[3][16] = -0.532709; w1[3][17] = 1.244082; w1[3][18] = -0.264233;
		w1[4][0] = 2.263374; w1[4][1] = 1.023849; w1[4][2] = -0.331190; w1[4][3] = -1.113271; w1[4][4] = -0.481672; w1[4][5] = 0.612255; w1[4][6] = 0.885745; w1[4][7] = 0.277659; w1[4][8] = 0.152636; w1[4][9] = -0.759591; w1[4][10] = -1.036185; w1[4][11] = 1.417303; w1[4][12] = 1.399967; w1[4][13] = 0.455211; w1[4][14] = -0.276722; w1[4][15] = 0.223897; w1[4][16] = 0.546051; w1[4][17] = -2.160171; w1[4][18] = -1.299445;
		w1[5][0] = 0.518210; w1[5][1] = 1.810684; w1[5][2] = 1.717064; w1[5][3] = -1.109350; w1[5][4] = 1.808813; w1[5][5] = 1.497307; w1[5][6] = -0.164426; w1[5][7] = -1.249963; w1[5][8] = 1.594136; w1[5][9] = -0.484002; w1[5][10] = 1.172055; w1[5][11] = 0.438501; w1[5][12] = 1.077568; w1[5][13] = -0.636005; w1[5][14] = -0.445795; w1[5][15] = -0.927332; w1[5][16] = -1.486168; w1[5][17] = 1.167487; w1[5][18] = 0.746513;
		w1[6][0] = -0.468745; w1[6][1] = 1.043974; w1[6][2] = -1.640883; w1[6][3] = 0.490878; w1[6][4] = -0.045497; w1[6][5] = -1.793901; w1[6][6] = -0.661866; w1[6][7] = 0.871435; w1[6][8] = 1.354049; w1[6][9] = 0.216862; w1[6][10] = 0.083540; w1[6][11] = 1.415415; w1[6][12] = -1.463577; w1[6][13] = -0.203821; w1[6][14] = -0.411375; w1[6][15] = -0.620165; w1[6][16] = 0.245724; w1[6][17] = 1.760054; w1[6][18] = -0.836105;


		w2[0] = -0.052509; w2[1] = -0.805185; w2[2] = -1.762332; w2[3] = 0.368993; w2[4] = -2.860540; w2[5] = -2.847379; w2[6] = 1.722926; w2[7] = -0.129739; 

		float *score = (float *)malloc(nfram*sizeof(float));
		for (int n=0;n<nfram;n++)
		{
			float hid[7];
			for (int i=0;i<7;i++) {hid[i]=0;}
			for (int i=0;i<7;i++)
			{
				for (int j=0;j<19;j++)
				{
					if (0==j) hid[i] += w1[i][j];
					else hid[i] += specGram[n][j-1]*w1[i][j];
				}
				hid[i] = 1/(1+exp(-hid[i]));
			}
			float out=0;
			for (int i=0;i<8;i++)
			{
				if(0==i) out += w2[i];
				else out += hid[i-1]*w2[i];
			}
			out = 1/(1+exp(-out));
			score[n] = out;
		}
		/*FILE *fpout = fopen("prob.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpout,"%f\n",score[n]);
		}
		fclose(fpout);*/
		//select candidata points
		float* output_norm_select = (float*)malloc(nfram * sizeof(float));
		//pick out candidate points
		for (int i=0;i<nfram;i++)
		{
			if (score[i]> THRESHOLD_8000)
			{
				output_norm_select[i] = score[i];
			}
			else
				output_norm_select[i] = 0;
		}

		/*FILE *fpSele = fopen("candiSele.txt","w");
		for (int n=0;n<nfram;n++)
		{
			fprintf(fpSele,"%f\n",output_norm_select[n]);
		}
		fclose(fpSele);*/

		//合并相邻的帧
		int flag = 0;
		float localMax = 0;
		int localPos = 0;
		for (int i=0;i<nfram;i++)
		{
			if (output_norm_select[i]>0)
			{
				if (flag>0)
				{
					if (output_norm_select[i]>localMax)
					{
						output_norm_select[localPos] = 0;
						localPos = i;
						localMax = output_norm_select[i];
					}
					else
					{
						output_norm_select[i] = 0;
					}
				}
				else
				{
					flag = 1;
					localPos = i;
					localMax = output_norm_select[i];
				}
			}
			else
			{
				flag = 0;
			}
		}
		//obtain localMaxPos
		int num=0;
		for (int i=0;i<nfram;i++)
		{
			if(output_norm_select[i]>0)
				num++;
		}
		if (num>0)
		{
			float *leapPositions = (float*)malloc(2*num * sizeof(float));
			num = 0;
			for (int i=0;i<nfram;i++)
			{
				if (output_norm_select[i]>0)
				{
					leapPositions[2*num] = i*0.01;
					leapPositions[2*num+1] = output_norm_select[i];
					num++;
				}
			}
			*len = num;
			return leapPositions;
		}
		else
		{
			*len = num;
			float *leapPositions = NULL;
			return leapPositions;
		}
	}


	//delete output_norm_select;

}