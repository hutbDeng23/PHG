/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: perf.c,v 1.57 2021/04/08 06:34:19 zlb Exp $
 *
 * Functions providing performance information */

#include "phg.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>	/* for strerror_r() */
#include <errno.h>

/* for open/read/close */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <signal.h>
#include <sys/time.h>
#ifdef __WIN__
# if HAVE_PSAPI
#  define FLOAT WIN_FLOAT
#  define CHAR WIN_CHAR
#  define BOOLEAN WIN_BOOLEAN
#  define COORD WIN_COORD
#  define VT_INT WIN_VT_INT
   /* Note: the FunctionEntry macro is used in psapi.h */
#  if defined(FunctionEntry)
#   undef FunctionEntry
#  endif
#  if defined(FunctionEntry_head)
#   undef FunctionEntry_head
#  endif
#  if defined(FunctionEntry_tail)
#   undef FunctionEntry_tail
#  endif
#  include <windows.h>
#  include <psapi.h>
#  undef FLOAT
#  undef CHAR
#  undef BOOLEAN
#  undef COORD
#  undef VT_INT
# endif	/* HAVE_PSAPI */
#elif defined(USETIMES)
# include <sys/times.h>
# include <unistd.h>
#else
# include <sys/resource.h>
#endif

/*---------------------------------------------------------------------------*/
/* Memory usage statistics */

/* Real time */
#define TimerSignal	SIGALRM
#define TimerType	ITIMER_REAL

/* Virtual time */
/*#define TimerSignal	SIGVTALRM
#define TimerType	ITIMER_VIRTUAL*/

/* Profiling time */
/*#define TimerSignal	SIGPROF
#define TimerType	ITIMER_PROF*/

static size_t	mem_peak = 0;
static size_t	mem_current = 0;
/* The sampling interval of the timer in ms, <= 0 disables the timer. */
static INT	sampling_interval = 200;
static char	statm[40];
static BOOLEAN	have_statm = FALSE;

static void
sample_mem_usage(int dummy)
{
    static BOOLEAN locked = FALSE;

#ifndef __WIN__

    int fd;
    struct rusage RU;
    size_t mem;

    if (locked)
	return;
    locked = TRUE;

    /* TODO: configure checks for different systems */

    /* first, try /proc/pid/statm for Linux systems */
    if (have_statm && (fd = open(statm, O_RDONLY)) != -1) {
	int rsspages;
	static char buffer[256];
	char *p = buffer;
	/* see linux-2.6.15/Documentation/filesystems/proc.txt */
	rsspages = read(fd, buffer, sizeof(buffer) - 1);
	close(fd);
	buffer[rsspages] = '\0';
#if 0
	/* Note: the following may hang on Fedora-16-i386 */
	sscanf(buffer, "%*d %d", &rsspages);
#else
	strtol(p, &p, 10);		/* first number () */
	rsspages = strtol(p, &p, 10);	/* second number */
#endif
	mem = (size_t)rsspages * (size_t)getpagesize();
    }
    else {
	/* next, try getrusage() */
	getrusage(RUSAGE_SELF, &RU);
#ifdef __APPLE_CC__
	mem = RU.ru_maxrss;
#else	/* defined(__APPLE_CC__) */
	mem = RU.ru_maxrss * (size_t)1024;
#endif	/* defined(__APPLE_CC__) */
	/*mem = RU.ru_maxrss * getpagesize();*/
    }

    if (mem > mem_peak)
	mem_peak = mem;
    mem_current = mem;

#else	/* !defined(__WIN__) */

#if HAVE_PSAPI
    HANDLE hProcess;
    PROCESS_MEMORY_COUNTERS pmc;

    if (locked)
	return;
    locked = TRUE;

    hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
			    FALSE, getpid());
    assert(hProcess != NULL);
    if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
#if 0
	printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
	printf( "\tPeakWorkingSetSize: 0x%08X\n",
		pmc.PeakWorkingSetSize );
	printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
	printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n",
		pmc.QuotaPeakPagedPoolUsage );
	printf( "\tQuotaPagedPoolUsage: 0x%08X\n",
		pmc.QuotaPagedPoolUsage );
	printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n",
		pmc.QuotaPeakNonPagedPoolUsage );
	printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n",
		pmc.QuotaNonPagedPoolUsage );
	printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage );
	printf( "\tPeakPagefileUsage: 0x%08X\n",
		pmc.PeakPagefileUsage );
#endif
	mem_current = pmc.WorkingSetSize;
	mem_peak = pmc.PeakWorkingSetSize;
    }
    CloseHandle( hProcess );
#endif	/* HAVE_PSAPI */
    
#endif	/* !defined(__WIN__) */

    locked = FALSE;
}

void
phgMemoryUsageReset(void)
{
    mem_peak = 0;
}

size_t
phgMemoryUsage1(GRID *g,
	size_t *current_min, size_t *current_avg, size_t *current_max,
	size_t *peak_min,    size_t *peak_avg,    size_t *peak_max)
{
    sample_mem_usage(0);

#if USE_MPI
    if (g != NULL && g->nprocs > 1 && !phgMasterSlave) {
	double a[2][3], b[2][3];
	a[0][0] = a[0][1] = a[0][2] = mem_current;
	a[1][0] = a[1][1] = a[1][2] = mem_peak;
	MPI_Allreduce(a, b, 2, MPI_3DOUBLE, MPI_MSM, g->comm);
	b[0][1] /= g->nprocs;
	b[1][1] /= g->nprocs;

	if (current_min != NULL)
	    *current_min = Roundsize_t(b[0][0]);
	if (current_avg != NULL)
	    *current_avg = Roundsize_t(b[0][1]);
	if (current_max != NULL)
	    *current_max = Roundsize_t(b[0][2]);

	if (peak_min != NULL)
	    *peak_min = Roundsize_t(b[1][0]);
	if (peak_avg != NULL)
	    *peak_avg = Roundsize_t(b[1][1]);
	if (peak_max != NULL)
	    *peak_max = Roundsize_t(b[1][2]);

	return (size_t)b[0][2];		/* return max(mem_current) */
    }
#else	/* USE_MPI */
    Unused(g);
#endif	/* USE_MPI */

    if (current_min != NULL)
	*current_min = mem_current;
    if (current_avg != NULL)
	*current_avg = mem_current;
    if (current_max != NULL)
	*current_max = mem_current;

    if (peak_min != NULL)
	*peak_min = mem_peak;
    if (peak_avg != NULL)
	*peak_avg = mem_peak;
    if (peak_max != NULL)
	*peak_max = mem_peak;


    return mem_current;
}

size_t
phgMemoryUsage(GRID *g, size_t *peak)
{
#if 0
    static BOOLEAN warn_omp = TRUE;

    if (phgMaxThreads > 1 && warn_omp) {
	warn_omp = FALSE;
	phgPrintf("WARNING: phgMemoryUsage is not thread-aware!\n");
    }
#endif
    return phgMemoryUsage1(g, NULL, NULL, NULL, NULL, NULL, peak);
}

size_t
phgMemoryPeakReset(size_t value)
{
    size_t old_value = mem_peak;
    mem_peak = value;

    return old_value;
}

size_t
phgMemoryPeakRestore(size_t value)
{
    size_t old_value = mem_peak;

    if (mem_peak < value)
	mem_peak = value;

    return old_value;
}

static BOOLEAN timer_initialized = FALSE;
#ifndef __WIN__
static struct itimerval new, old;
static struct sigaction sa, old_sa;
#endif	/* !defined(__WIN__) */

static void
timer_init(void)
{
    static BOOLEAN first_call = TRUE;
    int fd;

    if (first_call) {
	phgOptionsRegisterInt("-mem_sampling_rate", "Memory usage sampling "
		"rate (milliseconds, 0 disables memory sampling)",
		&sampling_interval);
	first_call = FALSE;
	return;
    }

    if (timer_initialized)
	return;
    timer_initialized = TRUE;

    sprintf(statm, "/proc/%d/statm", (int)getpid());
    if ((fd = open(statm, O_RDONLY)) != -1) {
	have_statm = TRUE;
	close(fd);
    }

    if (sampling_interval <= 0)
	return;

#ifndef __WIN__
    /* setup timer to sample memory usage */
    sa.sa_handler = sample_mem_usage;
    sigemptyset (&sa.sa_mask);
    /*sigfillset (&sa.sa_mask);*/
    sa.sa_flags = SA_RESTART;
    if (sigaction(TimerSignal, &sa, &old_sa))
	if (phgVerbosity > 0)
	    perror("sigaction 0");
    new.it_value.tv_sec = sampling_interval / 1000;
    new.it_value.tv_usec = (sampling_interval % 1000) * 1000;
    new.it_interval = new.it_value;
    if (setitimer(TimerType, &new, &old))
	if (phgVerbosity > 0)
	    perror("setitimer 0");
#endif	/* !defined(__WIN__) */

    return;
}

static void
timer_finalize(void)
{
    if (!timer_initialized)
	return;
    timer_initialized = FALSE;
 
#ifndef __WIN__
    /* restore timer */
    if (setitimer(ITIMER_REAL, &old, NULL))
	if (phgVerbosity > 0)
	    perror("setitimer 1");
    if (sigaction(SIGALRM, &old_sa, NULL))
	if (phgVerbosity > 0)
	    perror("sigaction 1");
#endif	/* !defined(__WIN__) */

    return;
}

/*---------------------------------------------------------------------------*/
/* Floating point operations */

static BOOLEAN use_perfctr = TRUE;

#if USE_PAPI

#include <papiStdEventDefs.h>
#include <papi.h>

#if USE_OMP
# include <pthread.h>		/* for pthread_self() */
#endif	/* USE_OMP */

static int (*func)(float *rtime, float *ptime, long_long *count, float *mf)
	= NULL;

static long_long count = 0;
static float rtime, ptime, mf;

static double time0 = -1., elapsed_time, last_time;
static long_long last_count = 0;

static void
set_time_count(void)
{
    double tarray[3];
    long_long tmp;

    phgGetTime(tarray);
    if (time0 < 0)
	time0 = tarray[2];
    elapsed_time = tarray[2] - last_time;
    tmp = count - last_count;
    last_count = count;
    last_time = tarray[2];
    count = tmp;
}

static void
papi_error(const char *funcname, int code)
{
#if !defined(PAPI_VERSION_MAJOR) || PAPI_VERSION_MAJOR(PAPI_VERSION) < 5
    char errmsg[PAPI_MAX_STR_LEN], *sep = ": ";
    PAPI_perror(code, errmsg, sizeof(errmsg));
#else	/* PAPI_VERSION_MAJOR */
    char *errmsg = "", *sep="";
    Unused(code);
    PAPI_perror(NULL);
#endif	/* PAPI_VERSION_MAJOR */
#if 0
    fprintf(stderr, "Initialization of PAPI_flops failed: %s", errmsg);
    if (code == PAPI_ESYS)
	perror(" ");
    else
	fprintf(stderr, "\n");
#else
    /*if (code == PAPI_ESYS) {
	char errmsg1[PAPI_MAX_STR_LEN];
	phgInfo(0, "Initialization of %s failed: %s %s\n", funcname, errmsg,
		strerror_r(errno, errmsg1, sizeof(errmsg1)));
    }
    else*/ {
	phgInfo(0, "Initialization of %s failed%s%s.\n", funcname, sep, errmsg);
    }
#endif
}

#ifdef EVENT
# undef EVENT
#endif
#if (SIZEOF_PHG_FLOAT == SIZEOF_DOUBLE)
# define EVENT		PAPI_DP_OPS
# define EVENT_NAME	"PAPI_DP_OPS"	/* PAPI_VEC_DP */
#elif (SIZEOF_PHG_FLOAT == SIZEOF_FLOAT)
# define EVENT		PAPI_SP_OPS
# define EVENT_NAME	"PAPI_SP_OPS"	/* PAPI_VEC_SP */
#endif

#if 0		/*---------------------------------------------------*/
/* This is a workaround for FPE in PAPI_get_{real,virt}_usec on LSSC3 */
static long long
get_time(int which)
{
    double t[3];
    phgGetTime(t);
    return (long long)(t[which] * 1e6);
}
#define PAPI_get_real_usec()	get_time(2)
#define PAPI_get_virt_usec()	get_time(0)
#endif		/*---------------------------------------------------*/

#ifdef EVENT
static int
PAPI_vflops(float *rtime, float *ptime, long long *flpops, float *mflops)
{
    static int initialized = 0;
    static long long total = 0, initial_rt = 0, initial_pt = 0, last_pt = 0;
    long long rt, pt;
    int retval;
    int Events[] = {EVENT};
    long long value = 0;

    if (!initialized) {
	initialized = 1;
	initial_rt = PAPI_get_real_usec();
	last_pt = initial_pt = PAPI_get_virt_usec();
	total = 0;

	return PAPI_start_counters(Events, 1);
    }

    if ((retval = PAPI_stop_counters(&value, 1)) != PAPI_OK) {
	initialized = 0;
	return retval;
    }

    rt = PAPI_get_real_usec();
    pt = PAPI_get_virt_usec();

    *flpops = (total += value);

    *rtime = (rt - initial_rt) * 0.000001;
    *ptime = pt - initial_pt;
    last_pt = pt - last_pt;
    *mflops = (last_pt != 0.0 ? value / last_pt : 0.0);
    *ptime *= 0.000001;

    last_pt = pt;

    if ((retval = PAPI_start_counters(Events, 1)) != PAPI_OK)
	initialized = 0;
    return retval;
}
#endif	/* defined(EVENT) */

static void
papi_init(void)
{
    int code;
    static int flag = -1;

    if (flag < 0) {
	/* register options */
	phgOptionsRegisterNoArg("-use_perfctr", "Use performance counters",
				&use_perfctr);
	flag = 0;
	return;
    }

    if (flag > 0 || !use_perfctr)
	return;

    flag = 1;

    code = PAPI_library_init(PAPI_VER_CURRENT);
    if (code != PAPI_VER_CURRENT) {
	if (phgVerbosity > 0)
	    papi_error("PAPI_library_init", code);
	return;
    }

#if USE_OMP
    code = PAPI_thread_init(pthread_self);
    if (code != PAPI_OK) {
	if (phgVerbosity > 0)
	    papi_error("PAPI_thread_init", code);
	return;
    }
#endif	/* USE_OMP */

#if 0
    PAPI_set_domain(PAPI_DOM_ALL);
    PAPI_set_granularity(/*PAPI_GRN_SYS*/PAPI_GRN_SYS_CPU);
#endif

#if 0
    func = PAPI_flops;
    code = func(&rtime, &ptime, &count, &mf);
    if (code != PAPI_OK) {
	if (phgVerbosity > 0)
	    papi_error("PAPI_flops", code);
	/*fprintf(stderr, "Trying PAPI_flips ...\n");*/
	func = PAPI_flips;
	code = func(&rtime, &ptime, &count, &mf);
	if (code != PAPI_OK) {
	    if (phgVerbosity > 0)
		papi_error("PAPI_flips", code);
	    func = NULL;
	}
    }
    if (phgVerbosity > 0) {
	if (func != NULL)
	    phgPrintf("PAPI: using %s.\n",
			func == PAPI_flops ? "PAPI_flops": "PAPI_flips");
	else
	    phgPrintf("PAPI: initialization failed.\n");
    }
#else
    if (PAPI_query_event(PAPI_FP_INS) == PAPI_OK) {
	if (phgVerbosity > 0)
	    phgPrintf("PAPI: using PAPI_FP_INS.\n");
	func = PAPI_flips;
    }
    else if (PAPI_query_event(PAPI_FP_OPS) == PAPI_OK) {
	if (phgVerbosity > 0)
	    phgPrintf("PAPI: using PAPI_FP_OPS.\n");
	func = PAPI_flops;
    }
#ifdef EVENT
    else if (PAPI_query_event(EVENT) == PAPI_OK) {
	if (phgVerbosity > 0)
	    phgPrintf("PAPI: using %s.\n", EVENT_NAME);
	func = PAPI_vflops;
    }
#endif	/* defined(EVENT) */
    else {
	if (phgVerbosity > 0)
	    phgPrintf("PAPI: FP counters unavailable on this platform.\n");
	func = NULL;
    }
    PAPI_shutdown();
#endif

    if (func != NULL) {
	set_time_count();
#if USE_OMP
	/* initialize PAPI's internal counters for all threads */
	{
	    long_long cnt;
#pragma omp parallel private(cnt, rtime, ptime, mf)
	    func(&rtime, &ptime, &cnt, &mf);
	}
#endif	/* USE_OMP */
    }
}

static void
papi_finalize(void)
{
    return;
}

double
phgPerfGetMflops(GRID *g, double *aggregate, double *since_last_call)
/* 'aggregate' = Mflops since the start of the program,
 * 'since_last_call' = Mflops since last call to this function.
 *
 * returns '*since_last_call' */
{
    double mflops[2];

    if (func == NULL) {
	if (aggregate != NULL)
	    *aggregate = 0.;
	if (since_last_call != NULL)
	    *since_last_call = 0.;
	return 0.;
    }

    if (phgMaxThreads == 1) {
	func(&rtime, &ptime, &count, &mf);
    }
#if USE_OMP
    else {
	long_long cnt = 0;
	count = 0;
#pragma omp parallel default(shared) \
	private(cnt, rtime, ptime, mf) reduction(+:count)
	{
	    func(&rtime, &ptime, &cnt, &mf);
	    /* printf("thread %d: cnt = %lld\n", phgThreadId, cnt);*/
	    count += cnt;
	}
    }
#endif	/* USE_OMP */

    set_time_count();

    mflops[0] = (elapsed_time == 0.0 ? 0. : count / elapsed_time * 1e-6);
    elapsed_time = last_time - time0;
    count = last_count;
    mflops[1] = (elapsed_time == 0.0 ? 0. : count / elapsed_time * 1e-6);

#if USE_MPI
    if (g != NULL && g->nprocs > 1 && !phgMasterSlave) {
	double a[2];
	MPI_Allreduce(mflops, a, 2, MPI_DOUBLE, MPI_SUM, g->comm);
	mflops[0] = a[0];
	mflops[1] = a[1];
    }
#endif	/* USE_MPI */

    if (aggregate != NULL)
	*aggregate = mflops[0];

    if (since_last_call != NULL)
	*since_last_call = mflops[1];

    return mflops[0];
}

#else	/* USE_PAPI */

static void
papi_init(void)
{
    static int flag = -1;

    if (flag < 0) {
	/* register options */
	phgOptionsRegisterNoArg("-use_perfctr", "Use performance counters",
				&use_perfctr);
	flag = 0;
	return;
    }
}

static void
papi_finalize(void)
{
    return;
}

double
phgPerfGetMflops(GRID *g, double *aggregate, double *since_last_call)
{
    if (aggregate != NULL)
	*aggregate = 0.;
    if (since_last_call != NULL)
	*since_last_call = 0.;

    return 0.;
}

#endif	/* USE_PAPI */

/*---------------------------------------------------------------------------*/
/* Processor speed */

#if HAVE_REGEX_H
/* regex(). Note: for MinGW requires mingw32-libgnurx or mingw64-libgnurx */
#include <string.h>
#include <sys/types.h>
#include <regex.h>

/* isspace() */
#include <ctype.h>
#endif	/* HAVE_REGEX_H */

static char *ps_file = NULL;
static char *ps_str = NULL;

static FLOAT proc_speed = 1., ps_min = 1., ps_avg = 1., ps_max = 1.;
static BOOLEAN ps_same = TRUE;

static void
ps_init(void)
/* Note: the format of the processor speed specifications can be either:
 *	regular expression: speed
 *	regular expression: speed
 *	......
 * (new-line separated), or:
 * 	regular expression: speed; regular expression: speed; ... ...
 * (semicolon separated). */
{
    static BOOLEAN first_call = TRUE;

#if USE_MPI
    int ret;
    FILE *f;
    char *str = NULL;
    double msm[6];
#if HAVE_REGEX_H
    regex_t re;
    double speed;
    char *p, *q, pattern[128], pname[MPI_MAX_PROCESSOR_NAME] = "";
#endif	/* HAVE_REGEX_H */
#endif	/* USE_MPI */

    if (first_call) {
	first_call = FALSE;
	phgOptionsRegisterFilename("-proc_speed_file", "Processor speed file",
					&ps_file);
	phgOptionsRegisterString("-proc_speed_spec", "Processor speed specs",
					&ps_str);
	phgOptionsRegisterFloat("-proc_speed", "Processor speed", &proc_speed);
	return;
    }

#if USE_MPI
    if (phgNProcs == 1)
	return;

    if (ps_file != NULL && (f = fopen(ps_file, "r")) != NULL) {
	/* check process speed file first. */
	size_t size;
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	str = phgAlloc(size + 1);
	fseek(f, 0, SEEK_SET);
	if (size > 0 && fread(str, size, 1, f) != 1)
		size = 0;
	str[size] = '\0';
	fclose(f);
    }
    else {
	/* check process speed string next. */
	str = ps_str;
    }

    if (str == NULL)
	goto post;

#if HAVE_REGEX_H

    for (p = str; *p != '\0'; p++)
	    if (*p == '\n')
		*p = ';';
    p = str;
	
    while (TRUE) {
	while (isspace(*p))
	    p++;
	if (*p == '\0')
	    break;
	/* read the pattern */
	q = strchr(p, ':');
	if (q == NULL)
	    phgError(1, "error in processor speed specifications.\n");
	assert(q - p < sizeof(pattern));
	*q = ' ';
	sscanf(p, "%s", pattern);
	/* read the speed */
	p = q + 1;
	q = strchr(p, ';');
	if (q == NULL)
	    q = p + strlen(p);
	else
	    *q = ' ';
	sscanf(p, "%lf", &speed);
	p = q;
	/* match processor name against the pattern */
	ret = regcomp(&re, pattern, 0 /*| REG_EXTENDED | REG_ICASE*/);
	if (ret != 0) {
	    char msg[128];
	    regerror(ret, &re, msg, sizeof(msg));
	    phgWarning("Error in the pattern \"%s\": %s\n", pattern, msg);
	    continue;
	}
	if (pname[0] == '\0')
	    MPI_Get_processor_name(pname, &ret);
	ret = regexec(&re, pname, (size_t)0, NULL, 0);
	regfree(&re);
#if 0
	phgInfo(-1, "pattern=\"%s\", speed=%lf %c\n", pattern, speed,
			ret == 0 ? '*' : ' ');
#endif
	if (ret == 0) {
	    proc_speed = speed;
	    break;
	}
    }

#else	/* HAVE_REGEX_H */

    phgError(1, "can't parse processor speed specs (regex() unavailable).\n");

#endif	/* HAVE_REGEX_H */

    if (str != ps_str)
	phgFree(str);

post:
    assert(proc_speed > 0);

    /* compute min/avg/max speed */
    msm[0] = msm[1] = msm[2] = proc_speed;
    MPI_Allreduce(msm, msm + 3, 1, MPI_3DOUBLE, MPI_MSM, phgComm);
    ps_min = msm[3];
    ps_avg = msm[4] / phgNProcs;
    ps_max = msm[5];

    ret = ((ps_max - ps_min) / ps_max <= 1e-3 ? TRUE : FALSE);
    MPI_Bcast(&ret, 1, MPI_INT, 0, phgComm);
    ps_same = ret;

    if (phgVerbosity > 0) {
	if (ps_same)
	    phgPrintf("*** All processors have the same speed %lg\n",
			(double)proc_speed);
	else
	    phgPrintf("*** Processors speed: min = %lf, avg = %lf, max = %lf\n",
			(double)ps_min, (double)ps_avg, (double)ps_max);
    }
#endif	/* USE_MPI */
}

FLOAT
phgGetProcessorSpeed(FLOAT mam[])
{
    if (mam != NULL) {
	mam[0] = ps_min;
	mam[1] = ps_avg;
	mam[2] = ps_max;
    }

    return proc_speed;
}

BOOLEAN
phgSameProcessorSpeed(void)
{
    return ps_same;
}

/*---------------------------------------------------------------------------*/

void
phgPerfInit(void)
{
    timer_init();
    papi_init();
    ps_init();
}

void
phgPerfFinalize(void)
{
    timer_finalize();
    papi_finalize();
}
