#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>

#ifndef BUFSIZE
#define BUFSIZE 	512
#endif

#ifndef PERM
#define PERM      0644
#endif

#ifndef	EPSON
#define	EPSON		0.01
#endif


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Reading one line from file
	return value indicate the characteristic of the read line
	isread = 0: end of file
	isread = 1: normal line
	isread = 2: blank line
	isread = 3: % line                                                         */
/*----------------------------------------------------------------------------*/
int readline(FILE* fp, char* oneline)
{
	int lcnter;
	int i;
	char header;
	int isread;
	int isblankl;
	
	// Read the whole line into "oneline" if only not reach file end.
	lcnter = 0;
	isread = 0;
	while(1) {
		
		if (feof(fp)) return isread;
		
		header = fgetc(fp);
		if (header != 10 && header != 13) {  /*-- 10=line feed, 13=carriage return --*/
			oneline[lcnter] = header;
			lcnter++;
		}
		else {
			oneline[lcnter] = '\0';
			break;
		}	
	}	
		
	// Check the "oneline" so as to determine reture value	
	// Neglect all beginning white spaces
	// isread = 0: end of file
	// isread = 1: normal line
	// isread = 2: blank line
	// isread = 3: comment line
	
	//printf("%s\n", oneline);
	if (strlen(oneline) == 0) return isread = 2;
	
	isblankl = 1;
	for (i = 0; i < (int) strlen(oneline); i++) {
		if (oneline[i] != 32 && oneline[i] != 9) {
			isblankl = 0;
			break;
		}
	}
	if (isblankl) return isread = 2;
			
	for (i = 0; i < (int) strlen(oneline); i++) {
		if (oneline[i] == 10 || oneline[i] == 11 || oneline[i] == 12 || oneline[i] == 13) {
			isread = 2;
			break;
		}
		else if (oneline[i] == '%' || oneline[i] == '#' || oneline[i] == '*') {
			isread = 3;
			break;
		}
		else { 
			isread = 1;
			break;
		}
	}
	
	return isread;		
}



/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* search "inchar" from end until reach "seq"_th "token"
	"inchar" are divided into two segments
	"charp1" from beginning to token (including token)
	"charp2" from token+1 to end (not including token)                         */
/*----------------------------------------------------------------------------*/
int segchar(char* inchar, char* charp1, char* charp2, char token, int seq)
{
	int isfound = 0;
	int i;
	int idir = 0;
	
	for (i = strlen(inchar) - 1; i >= 0; i--) {
		if (inchar[i] == token) {
			isfound++;
		}
		if (isfound == seq)	{
			idir = i+1;
			break;
		}
	}
	
	if (!isfound) {
		printf("\nNull Output from find_comm.\n");
		return 0;
	}
	
	for (i = 0; i < idir; i++) 
		charp1[i] = inchar[i];
	charp1[idir] = '\0';
	
	for (i = idir; i < (int) strlen(inchar); i++) 
		charp2[i-idir] = inchar[i];
	charp2[strlen(inchar)-idir] = '\0';

	return 1;		  
}



/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* copy file contents in name1 to name2                                       */
/*----------------------------------------------------------------------------*/
int copyfile(const char *name1, const char *name2)
{
	int infile, outfile;
	ssize_t nread;
	char buffer[BUFSIZE];
	
	if ((infile = open(name1, O_RDONLY)) == -1) return (-1);
	
	if ((outfile = open(name2, O_WRONLY|O_CREAT|O_TRUNC, PERM)) == -1) {
		close(infile);
		return (-1);
	}
	
	/*----- Read from name1 BUFSIZE chars at a time ---------*/
	while((nread = read(infile, buffer, BUFSIZE)) > 0) {
		/*----- Write buffer to outfile ------*/
		if (write(outfile, buffer, nread) < nread) {
			close(infile);
			close(outfile);
			return (-3);
		}
	}
	
	close(infile);
	close(outfile);
	
	if (nread == -1)
		return (-4);         /*---- error on last read ---*/
	else
		return (0);		       /*---- all is well ----*/

}		



/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* rm all files under dir named "name"                                        */
/*----------------------------------------------------------------------------*/
int rmallfiles(const char *name)
{
	struct dirent *d;
	DIR *dp;
	char pathname[256];
	
/* open the dir and check for failure */
	if ((dp = opendir(name)) == NULL) return -1;
	
	while (d = readdir(dp)) {
		if (d->d_ino != 0 && strcmp(d->d_name, ".") && strcmp(d->d_name, "..") ) {
			sprintf(pathname, "%s/%s", name, d->d_name);
			if (remove(pathname) != 0) {
				fprintf(stderr, "Fail to rm %s\n\n", pathname);
				return -2;
			}
		}
	}				
					
	return 0;
}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* replace oldv with newv in data                                             */
/*----------------------------------------------------------------------------*/
int replace(char *data, const char *oldv, char *newv)
{
	int i, j;
	int lendata;
	int lenoldv;
	int lennewv;
	int lenndt;
	int tag = 0;
	int ntops;
	int isfound = 0;
	char* tmpdata;

	lendata = strlen(data);
	lenoldv = strlen(oldv);
	lennewv = strlen(newv);
	
	tmpdata = new char[lendata+lennewv];
	
	for (i = 0; i < lendata; i++) {
		for (j = 0; j < lenoldv; j++) {
			if (oldv[j] == data[i+j]) {
				tag = i+j;
				isfound = 1;
			}
			else {	
				isfound = 0;
				break;
			}
		}
		if (isfound) {
			break;
		}
	}
	
	ntops = tag - lenoldv + 1;
	strncpy(tmpdata, data, ntops);
	strcat(tmpdata, newv);
	lenndt = strlen(tmpdata);	
	//printf("tmpdata: %s\n", tmpdata);
	//printf("lentmpdata = %d\n", lenndt);				
	//printf("tag = %d\n", tag);				
	//printf("lendata = %d\n", lendata);				
	for (i = tag+1; i < lendata; i++) {
		//printf(": %c\n", data[i]);			
		//printf("lentmpdata = %d\n", lenndt+i-tag);	
		tmpdata[lenndt+i-tag-1] = data[i];
	}
	//printf("lentmpdata = %d\n", lendata+lennewv-lenoldv);	
	//printf("lentmpdata = %d\n", strlen(tmpdata));				
	tmpdata[lendata+lennewv-lenoldv] = '\0';
	
	strcpy(data, tmpdata);
	delete [] tmpdata;
	
	return 0;
}



/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Print error massage                                                        */
/*----------------------------------------------------------------------------*/
void _errmsg_(char *msg)
{
	if (msg == NULL) fprintf(stderr, "\nError: NULL string.\n\n");
	else fprintf(stderr, "\nError in processing '%s'.\n\n", msg);
	exit(-1);
}
	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Calculate the number of segments from minVal to maxVal in step
	maxVal may be changed if cannot evenly divided                             */
/*----------------------------------------------------------------------------*/
int set_n(float minVal, float &maxVal, float step)
{
	int n;
	float rbuf;

	if (step < 1.0e-16) {
		n = 1;
		maxVal = minVal;
	}	 
	else {
		rbuf = (maxVal - minVal) / step;
		if ((rbuf - (int)rbuf) > EPSON * step) {
			n = (int)ceil(rbuf) + 1;
		}
		else{
			n = (int)rbuf  + 1;
		}	
		maxVal = minVal + (n - 1) * step;
	}
	
	return(n);
}
		 	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Separate string "str" into string array "av" deliminated by "token"
	"av" points to different location of "str", and 
	"str" is destroyed                                                         */
/*----------------------------------------------------------------------------*/
int _ss2token_(char *str, char **av, char token)
{
	int ac = 0;
	int flag  = 0;
	char ptr;

	for (ptr=0; str[ptr]!= 0x0; ptr++) {
		if (str[ptr] == token) {
			str[ptr] = 0x0;
			if (flag == 1)   flag = 0;
		}
		else {
			if (flag == 0) {
				*(av+ac) = &str[ptr];
				ac++;
				flag = 1;
			}
		}
	}
	
	return(ac);
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Separate string "stro" into string array "av" deliminated by "token"
	"av" points to different location of "stro", and 
	Operate on a copy of "stro", so "stro" is not destroyed                    */
/*----------------------------------------------------------------------------*/
int _ssCtoken_(char *stro, char **av, char token)
{
	int ac = 0;
	int flag  = 0;
	char ptr;
	char *str;

	str = new char[strlen(stro)+1];
	strcpy(str, stro);

	for (ptr=0; str[ptr]!= 0x0; ptr++) {
		if (str[ptr] == token) {
			str[ptr] = 0x0;
			if (flag == 1)   flag = 0;
		}
		else {
			if (flag == 0) {
				*(av+ac) = &str[ptr];
				ac++;
				flag = 1;
			}
		}
	}
	
	//delete[] str;
	return(ac);
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Check the path from root and make dir along the path
	the path chars may be destroyed           			                        */
/*----------------------------------------------------------------------------*/
int mkpath(char *str, char token)
{
	int _ssCtoken_(char *str, char **av, char token);
	
	char *av[64];  
	char fodir[256];
	int ac = 0;
	int mdt;
	int i;
	DIR* dp;
	
	ac = _ssCtoken_(str, av, token);
	//for (i = 0; i < ac; i++) {
	//	printf("av[%d] = '%s', strlen = %d\n", i, av[i], strlen(av[i]));
	//}
	//av[ac-1][strlen(av[ac-1])] = '\0';	
	//printf("av[%d] = %s, strlen = %d\n", ac-1, av[ac-1], strlen(av[ac-1]));
	
	if (ac == 0) {
		if ((dp = opendir(av[0])) == NULL) {
			if ((mdt = mkdir(av[0], 0777)) != 0) {
				fprintf(stderr, "\nFail create dir %s code = %d!\n\n", av[0], mdt);
				return -1;
			}
		}
		else {
			closedir(dp);
			return 0;
		}				
	}
	
	sprintf(fodir, "/%s/", av[0]);
	for (i = 0; i < ac; i++) {
		if ((dp = opendir(fodir)) == NULL) {
			if ((mdt = mkdir(fodir, 0777)) != 0) {
				fprintf(stderr, "\nFail create dir %s code = %d!\n\n", fodir, mdt);
				return -1;
			}
		}
		else {
			closedir(dp);
		}							

		if (i < ac - 1) {
			strcat(fodir, av[i+1]);
			strcat(fodir, "/");
		}
		//printf("fodir = '%s', length = %d\n", fodir, strlen(fodir));	
	}	
		
  return(0);
}

						
/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* zero a float array                                                         */
/*----------------------------------------------------------------------------*/
int zero(float *x, int n)
{
	int i;

	for (i = 0; i < n; i++) x[i] = 0.0;

	return 0;
}


