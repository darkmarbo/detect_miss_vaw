#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define M  15
#define N  40
#define E  500000000

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
int SAMFREQ;
int ReadFile(char *wfile)
{
	bool oflag=false;
	FILE *fp=NULL;
	WAVEHEAD head;
	SAMFREQ=-1;
	int sample_count=0,channel_num=0,readflag=0,nbits=0;
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
			nbits = head.psbits;

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

float * missDetect(char * filename,int *numFrame);
int main(int argc, char *argv[])
{
	FILE* fp;
	if ((fp = fopen(argv[2], "a+")) == NULL)
	{
		printf("Please provide the second parameter...");
		return 0;
	}

	char* filename = argv[1];
	int freq = SAMFREQ;
	int missNumFrame;
	float* missPositions = missDetect(filename,&missNumFrame);//单位是frame
	fprintf(fp, "%s\t", filename);

	if (missNumFrame == 0)
	{
		fprintf(fp, "0\t\n");
		return 0;
	}
	fprintf(fp, "1\t");
	fprintf(fp, "%d\t", missNumFrame);

	for(int n=0; n<missNumFrame; n++)
	{
		//fprintf(fp, "miss number %d located at: %d min %f sec\n", 
			//n,int(missPositions[n]),60*(missPositions[n]-int(missPositions[n])));
		fprintf(fp, "%.3f\t", missPositions[n]*60);
	}

	fprintf(fp, "\n");
	fclose(fp);
	//printf("done!\n");
	return 0;

}

float * missDetect(char * filename,int *numFrameRet)
{
	int fLength;
	long sigLength = ReadFile(filename);
	if (sigLength<1)
	{
		return NULL;
	}

	int winLength = SAMFREQ*0.02;
	int numFrame = int(sigLength/winLength*2)-1;
	int *flagSaver = (int *)malloc(numFrame*sizeof(int));
	for (int i=0;i<numFrame;i++) flagSaver[i] = 0;
	for (int i=0;i<numFrame;i++)
	{
		int zero_counter = 0;
		float eng = 0;
		for (int j=i*(winLength/2);j<i*(winLength/2)+winLength-1;j++)
		{
			if(0==allbuf[j]) zero_counter++;
			eng += allbuf[j]*allbuf[j];
		}
		if(zero_counter>M && zero_counter<N && eng>E) 
			flagSaver[i]=1;
		zero_counter = 0;
		eng = 0;
	}
	int tic;
	for (int i=0;i<numFrame;i++)
	{
		tic = 0;
		if(flagSaver[i]>0) 
		{
			if (tic>0) flagSaver[i] = 0;
			else tic++;
		}
	}

	int missNumFrame = 0;
	for (int i=0;i<numFrame;i++)
	{
		if(flagSaver[i]>0) 
		{
			missNumFrame++;
		}
	}
	*numFrameRet = missNumFrame;
	float* effPositions = (float*)malloc(missNumFrame*sizeof(float));
	int tick = 0;
	for (int i=0;i<numFrame;i++)
	{
		if(flagSaver[i]>0) 
		{
			effPositions[tick++] = (float(i)/float(6000));
		}
	}
	return effPositions;
}