#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio>


/******************************************************************************/
/* Generate H-k and S-function file name with path                            
	opt = 1: Use old name format: hk2dx7.AXX, hk2dx7.AXX.00..., sfr2dx7.AXX
	opt = 2: Use new name format: hk2d.HL.AXX.x7, hk2d.HL.AXX.x7.00..., 
				sfr2d.HL.AXX.x7   
	Make sure the path 'fodir' exists. 
	In this case the path 'fodir' has been made in inputting path name         */
/******************************************************************************/
int create_hkfname(char *fpth_file, char *fodir, char *hkfile, 
						const int& xc_flag, const int& opt)
{ 
	char seg1[16];
	char seg2[8];

	extern int segchar(char* inchar, char* charp1, char* charp2, char token, int seq);

	if (opt == 1) {
		if (segchar(hkfile, seg1, seg2, '.', 1) != 1) return -3;
		seg1[strlen(seg1)-1] = '\0';
		sprintf(fpth_file, "%s%sx%d.%s", fodir, seg1, xc_flag, seg2);	
	}
	else if (opt == 2) {
		sprintf(fpth_file, "%s%s.x%d", fodir, hkfile, xc_flag);			
	}

	return 0;
}	


/******************************************************************************/
/* Generate H-k and S-function file name with path 
	This version is for moving BAZ                            
	opt = 1: Use old name format only for sfr2d.b010x7.AXX
	opt = 2: Use new name format only for sfr2d.HL.AXX.x7.b010   
	Make sure the path 'fodir' exists. 
	In this case the path 'fodir' has been made in inputting path name         */
/******************************************************************************/
int create_hkfnameb(char *fpth_file, char *fodir, char *hkfile, 
						const float& baz_lb, const int& xc_flag, const int& opt)
{ 
	char seg1[16];
	char seg2[8];

	extern int segchar(char* inchar, char* charp1, char* charp2, char token, int seq);

	if (opt == 1) {
		if (segchar(hkfile, seg1, seg2, '.', 1) != 1) return -3;
		seg1[strlen(seg1)-1] = '\0';
		sprintf(fpth_file, "%s%s.b%03dx%d.%s", fodir, seg1, (int)baz_lb, xc_flag, seg2);	
	}
	else if (opt == 2) {
		sprintf(fpth_file, "%s%s.x%d.b%03d", fodir, hkfile, xc_flag, (int)baz_lb);			
	}

	return 0;
}	

