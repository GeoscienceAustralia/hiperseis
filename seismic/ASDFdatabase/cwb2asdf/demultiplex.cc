/********************************************************
 * Description:
 *   Reads multiplexed waveform data in mseed format, 
 *   dumped from Antelope, and outputs a number of 
 *   aggregated mseed files grouped by network and 
 *   station
 *   
 *  References:
 *    libmseed : https://github.com/iris-edu/libmseed/tree/master
 *
 * CreationDate:   03/04/19
 * Developer:      rakib.hassan@ga.gov.au
 *
 * Revision History:
 *   LastUpdate:     03/04/19   RH
 *   LastUpdate:     dd/mm/yyyy  Who     Optional description
 ********************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <map>
#include <fstream>
#include <unistd.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
    #include <libmseed.h>
#ifdef __cplusplus
}
#endif

typedef map< string, map<string, MSTraceGroup *> > TraceCollection ;
static double timetol     = -1.0; /* Time tolerance for continuous traces */
static double sampratetol = -1.0; /* Sample rate tolerance for continuous traces */


void process_mem_usage(double& vm_usage, double& resident_set)
{
   /********************************************************
    * Taken from: 
    * https://gist.github.com/thirdwing/da4621eb163a886a03c5
    ********************************************************/
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

/********************************************************
* Reads multiplexed mseed files and outputs aggregated 
* traces grouped by Network and Station
********************************************************/
int main(int argc, char **argv)
{
    /********************************************************
    * Process arguments:
    ********************************************************/
    if(argc < 4)
    {
        cout << endl << "Usage: ./demultiplex input_mseed output_folder output_prefix" << endl << endl;
        exit(0);
    }

    string fname   = string(argv[1]);
    string ofolder = string(argv[2]);
    string oprefix = string(argv[3]);
    MSRecord *msr  = NULL;
    int retcode;
    int verbose    = 0;
    
    cout << "Processing file " << fname << endl;    

    TraceCollection tc;

    while ( (retcode = ms_readmsr (&msr, fname.c_str(), -1, NULL, NULL, 1, 1, verbose)) == MS_NOERROR )
    {
        if(tc[msr->network][msr->station] == NULL)
        {
            tc[msr->network][msr->station] = mst_initgroup (NULL);
        }
        else
        {
            /********************************************************
            * Add record to trace-group that lives inside a nested
            * map keyed by network-code and then by station-code
            ********************************************************/

            MSTraceGroup *tg = tc[msr->network][msr->station];
            mst_addmsrtogroup(tg, msr, 0, timetol, sampratetol);
        }
    }
    cout << "Read all data.." << endl;

    TraceCollection::iterator tciter;
    int netsta_count = 0;
    int trace_count = 0;
    for (tciter = tc.begin(); tciter != tc.end(); tciter++)
    {
        for (map<string, MSTraceGroup*>::iterator tgiter=tciter->second.begin();
             tgiter != tciter->second.end(); tgiter++)
        {
            MSTraceGroup *mstg = tgiter->second;

            /* Combine traces */
            mst_groupheal(mstg, timetol, sampratetol);

            char ofn[512] = {0};
            sprintf(ofn, "%s/%s.%s.%s.mseed", 
                    ofolder.c_str(), oprefix.c_str(),
                    tciter->first.c_str(), tgiter->first.c_str());
            mst_writemseedgroup(mstg, ofn, 1, 4096, DE_STEIM2, 1, 1);

            netsta_count += 1;
            trace_count += mstg->numtraces;
        }
    }

    cout << "Wrote " << netsta_count << " mseed files, containing " << trace_count << " traces.." << endl;

    if ( retcode != MS_ENDOFFILE )
        ms_log (2, "Cannot read %s: %s\n", fname.c_str(), ms_errorstr(retcode));

    /* Cleanup memory and close file */
    ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);

    double vm, rss;
    process_mem_usage(vm, rss);
    cout << "Memory Used (MB) VM: " << vm/1024. << ", RSS: " << rss/1024. << endl;
}

