// probeDlg.cpp : implementation file
//

#include "stdafx.h"
#include "probe.h"
#include "probeDlg.h"
//*******************************************beg
#include <iostream>
#include <fstream>
#include "afxwin.h"		// библиотека MFC
#include "conio.h"
#include "math.h"
#include "time.h"
///////added
#include <stdio.h>
//#include <conio.h>
#include <stdlib.h>
#include <windows.h>
//#include <time.h>
#include <dos.h>
//#include <math.h>
#include <string.h>
///////\added
using namespace std;
//*******************************************end

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
//****************************************beg //definition//////////////////////////////////////////////////////////////////
int k,n,flc,csc,offset,hoffset_fazacf,show_corresp,show_num, interspike_prev,bcx,bcy,interspike[30000],binn,bin[100];
float interspike_sum, interspike_sum_max, i_float, csc_float, scale_rihist, scale_rihist_hor, show_red, rgb_scale;
int rihist_offset;
float bufsgm[20000][1200];
	double scale, offset_fazacf, scale_fazacf;
	double temp,sum,fazacf_prev,sum_prev;
	char strtemp[100];
	ifstream ifs;
	//////////added
	const int Nmax=1200; //1000;  //number of neurons
//const int tau=10; //in milliseconds
const int Mmax=500; //400; //(int)(100*10/tau); //number of samples
int M,N,wt,Nnb,Nnbc_rab;
int Nnbc, base, bin_integral;
const int Nnbmax=20000; //Buffer size
const int TavP=20;//7; //Buffer for averaging
const int TW=111111;//1010000;//A constant for eligibility phases
typedef float vect[Nmax]; //vectors definition
typedef int vecti[Nmax]; //integer vectors definition
int grphase[Nmax];// Phase of granule cells
int grocc[Nmax];  // ## of occurences of GrCs
int i_random;
vect swzi; // vectors  of connection weights
vecti vhodav;
vecti act; // vector of neuron activity; counter of occurence
vecti vhoda[Mmax];   //input sequence (adresa vozb. volokon)
vecti vhod[Mmax];   //input sequence  ("0" i "1")
float sgm[Nmax],oldsgm[Nmax],fx0sgm[Nmax],fx1sgm[Nmax],fx2sgm[Nmax]; //variables FOR connection weights
float hhi[100000];//connections changes function
float purk, purkin, btt, conI,intgr, sumwt,epsi, alf; //variables FOR
int Nss, i, s, q,tm, tm_show, tms, tminus, fazacf, fazat, Tpred, Tepsioff, Ipred, cht, Tmsp, fazan; //variables FOR
int cfexci, tmg, adrr,alfama,j,ik,uh,ic,ja,ib,iwa,iwe,iwb;
int elig_i, qq_check, logic_i, logic_failed;
int tm_ust1, tm_ust2;
float perem, purk_mult, intgr_mult;
void generation();// generation of a sequence of inputs
typedef char str[10];//string for filenames
str StAllt,StAt1,StSe1; //filenames
FILE *fall,*fall1,*fse1,*hall,*sall; //files for records
FILE *fnetww; //File of network weights
FILE *timcs_file, *sgm_file;
float bufs[Nnbmax]; //buffer dlia vyvoda
float buff[Nnbmax]; //buffer dlia vyvoda
int bufe[Nnbmax];//CF excitations
int buft[Nnbmax];//tm buffer
int FphaseVal,SphaseVal; // constant of changes in the first phase
								 // and constant of changes in the second phase
int Pfunc[5][Mmax];//values of function of phase
int fNo;//current # of the acting function
int Trefr, Texci, Refr;
int Ncs_pred, Ncs_prev;
int OutPu, PhaseP, limIP, TimpP;
float intpurk,LimTek;
int timcs_i;
char str1[300],str2[300],str3[300], str_temp[300];
int step;
float freq1, freq2;
int ntek,ttek[TavP];
int iaa,ibb,cons,TavPstoh;
float eligc[TW],eligc1[Nmax],eligc2[Nmax],eligc3[Nmax]; //Eligibility function; its 1 for finite number of t around 0, then 0
int timcs[10000000], timss[1000000], Ncs, ss, cs, tekss[1000],tekcs[10],bdip;//tek=tekushshee
short FqCS[10000000], Tf=1000;//Frequency(number) of CSs in Tf tacts of time
int TecFrCell;
float rab,rab1,rab10,rab11,rab12,rab20,rab21,rab22;
int irab,irab1,irab2,irab3;
int qq,a,iph=0;
int umka_max, noise_max;
int eligplus;
float fi_cylinder, r_cylinder;
double epsi_umenshitel; //float - up to 2^16 = 65536 (?), double - up to 2^32 = 4*10^9
void Demo();//Demo functsionirovaniya seti (na vyhod)
void StDemo();//Stationary version: Demo functsionirovaniya seti (na vyhod)
void Savo();//Save network connections
void Extc();//Extract connections for the network
void PrFqCS();//Printout of the CSs frequences
void other_PC_influence_calculation();
void renewing_of_weights_sampling_point (int tm_given);
void fprintf_and_printf();
void phase_control();
void unknown_function_1 ();
void initialization_of_everything();
void general_randomization(int rndmz);
void protecting_overflow_of_grphase();
void bufferization_of_purk_and_alf();
void show_gc_activity();
void show_gc_to_pc_weights();
void eligibility_as_grphase_initialization();
void generation2();
void elig_init_as_grphase();
void make_jump1();
void make_jump2();
void make_epsi_little();
void initialize_N_M_wt();
void initialize_wt();
void my_refresh_screen();
int main_calculation();
void str2_cfout_make();
void str2_timcs_make();
void str2_sgm_make();
void make_str2_cfout_fprintf();
void make_str2_timcs_fprintf();
void make_str2_sgm_fprintf();
float lplot_sum;
int gen; //codition to resume generation
int oneval;//a parameter in One Value eligibility
int sfc0,kfc0;//for control purposes
int Mcount_min;
int refresh_each_time;
int epsi_count, eligc_mode;
int razr;//Razriadka - constant for associations
int sinmode;
int Mcount, Mcount_max;
int Tequ, Tmin, Tmax, Tmin1, Tmin2; //Equilibrium Period
float OtherNC; //Amplitude of other PC influence on the IOC
int others; //equivalent action of others PCs
float OLDSQ, SQ, DST, rabs;
int iy,chis,Nav,Nav1;
int Tjump1, Tjump2, Ncount;
int allowexit1,allowexit2,allowfprintf;
int grn;
int epsi_umka, noise;
int progmode, create_file;
float epsi_basal, epsi_basal_first;
	//////////\added
	//************************************end
/////////////////////////////////////////////////////////////////////////////
// CProbeDlg dialog

CProbeDlg::CProbeDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CProbeDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CProbeDlg)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);

	//*****************************************beg // initialization///////////////////////////////////////////////////////
	Tpred=6000000;//900000000;
	n=10; offset=1; create_file=0;
	temp=sum=fazacf=fazacf_prev=csc=0.;
	scale=1.;
	offset_fazacf=0.; scale_fazacf=1.; hoffset_fazacf=0; scale_rihist = 1.; scale_rihist_hor = 1.; rihist_offset=0; show_red=0;
	show_corresp=0; show_num=0; refresh_each_time=0;
	//N=20; M=10; step=0; sinmode=0; noise=0; eligc_mode=2; epsi_count=2; epsi_umka=1; tm_show=1998000;
	N=300; M=200; step=0; sinmode=1; noise=0; eligc_mode=2; epsi_count=6; epsi_umka=1; tm_show=1998000; rgb_scale=1.;
	//for (int i=0; i<20000; i++) for (int j=0; j<12000; j++) bufsgm[i][j]=0.;
	//*****************************************end
}

void CProbeDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CProbeDlg)
	// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
	DDX_Control(pDX, IDC_SLIDER1, my_slider1);
}

BEGIN_MESSAGE_MAP(CProbeDlg, CDialog) 
	//{{AFX_MSG_MAP(CProbeDlg)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, OnButton1)
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BUTTON2, &CProbeDlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &CProbeDlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON4, &CProbeDlg::OnBnClickedButton4)
	ON_BN_CLICKED(IDC_BUTTON5, &CProbeDlg::OnBnClickedButton5)
	ON_BN_CLICKED(IDC_BUTTON6, &CProbeDlg::OnBnClickedButton6)
	ON_BN_CLICKED(IDC_BUTTON7, &CProbeDlg::OnBnClickedButton7)
	ON_BN_CLICKED(IDC_BUTTON8, &CProbeDlg::OnBnClickedButton8)
	ON_BN_CLICKED(IDC_BUTTON9, &CProbeDlg::OnBnClickedButton9)
	ON_BN_CLICKED(IDC_BUTTON10, &CProbeDlg::OnBnClickedButton10)
	ON_BN_CLICKED(IDC_BUTTON11, &CProbeDlg::OnBnClickedButton11)
	ON_BN_CLICKED(IDC_BUTTON12, &CProbeDlg::OnBnClickedButton12)
	ON_BN_CLICKED(IDC_BUTTON13, &CProbeDlg::OnBnClickedButton13)
	ON_BN_CLICKED(IDC_BUTTON14, &CProbeDlg::OnBnClickedButton14)
	ON_BN_CLICKED(IDC_BUTTON15, &CProbeDlg::OnBnClickedButton15)
	ON_BN_CLICKED(IDC_BUTTON16, &CProbeDlg::OnBnClickedButton16)
	ON_BN_CLICKED(IDC_BUTTON17, &CProbeDlg::OnBnClickedButton17)
	ON_BN_CLICKED(IDC_BUTTON18, &CProbeDlg::OnBnClickedButton18)
	ON_BN_CLICKED(IDC_BUTTON19, &CProbeDlg::OnBnClickedButton19)
	ON_BN_CLICKED(IDC_BUTTON20, &CProbeDlg::OnBnClickedButton20)
	ON_BN_CLICKED(IDC_BUTTON21, &CProbeDlg::OnBnClickedButton21)
	ON_BN_CLICKED(IDC_BUTTON22, &CProbeDlg::OnBnClickedButton22)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, &CProbeDlg::OnNMCustomdrawSlider1)
	ON_BN_CLICKED(IDC_BUTTON23, &CProbeDlg::OnBnClickedButton23)
	ON_BN_CLICKED(IDC_BUTTON24, &CProbeDlg::OnBnClickedButton24)
	ON_BN_CLICKED(IDC_BUTTON25, &CProbeDlg::OnBnClickedButton25)
	ON_BN_CLICKED(IDC_BUTTON26, &CProbeDlg::OnBnClickedButton26)
	ON_BN_CLICKED(IDC_BUTTON27, &CProbeDlg::OnBnClickedButton27)
	ON_BN_CLICKED(IDC_BUTTON28, &CProbeDlg::OnBnClickedButton28)
	ON_BN_CLICKED(IDC_BUTTON29, &CProbeDlg::OnBnClickedButton29)
	ON_BN_CLICKED(IDC_BUTTON30, &CProbeDlg::OnBnClickedButton30)
	ON_BN_CLICKED(IDC_BUTTON31, &CProbeDlg::OnBnClickedButton31)
	ON_BN_CLICKED(IDC_BUTTON32, &CProbeDlg::OnBnClickedButton32)
	ON_BN_CLICKED(IDC_BUTTON33, &CProbeDlg::OnBnClickedButton33)
	ON_BN_CLICKED(IDC_BUTTON34, &CProbeDlg::OnBnClickedButton34)
	ON_BN_CLICKED(IDC_BUTTON35, &CProbeDlg::OnBnClickedButton35)
	ON_BN_CLICKED(IDC_BUTTON36, &CProbeDlg::OnBnClickedButton36)
	ON_BN_CLICKED(IDC_BUTTON37, &CProbeDlg::OnBnClickedButton37)
	ON_BN_CLICKED(IDC_BUTTON38, &CProbeDlg::OnBnClickedButton38)
	ON_BN_CLICKED(IDC_BUTTON39, &CProbeDlg::OnBnClickedButton39)
	ON_BN_CLICKED(IDC_BUTTON40, &CProbeDlg::OnBnClickedButton40)
	ON_BN_CLICKED(IDC_BUTTON41, &CProbeDlg::OnBnClickedButton41)
	ON_BN_CLICKED(IDC_BUTTON42, &CProbeDlg::OnBnClickedButton42)
	ON_BN_CLICKED(IDC_BUTTON43, &CProbeDlg::OnBnClickedButton43)
	ON_BN_CLICKED(IDC_BUTTON44, &CProbeDlg::OnBnClickedButton44)
	ON_BN_CLICKED(IDC_BUTTON45, &CProbeDlg::OnBnClickedButton45)
	ON_BN_CLICKED(IDC_BUTTON46, &CProbeDlg::OnBnClickedButton46)
	ON_BN_CLICKED(IDC_BUTTON47, &CProbeDlg::OnBnClickedButton47)
	ON_BN_CLICKED(IDC_BUTTON48, &CProbeDlg::OnBnClickedButton48)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CProbeDlg message handlers
 
BOOL CProbeDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	// TODO: Add extra initialization here
	
	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CProbeDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
			my_slider1.SetRange(0,100,false);

		dc.MoveTo(100,300); dc.LineTo(1200,300);
		
	}
	else
	{
		CPaintDC dc(this); // device context for painting
				//********************************************************beg OnPaint////////////////////////////////////
		dc.MoveTo(100,300); dc.LineTo(1200,300);
		bcx=200; bcy=900;
		dc.MoveTo(bcx,bcy-100); dc.LineTo(bcx,bcy); dc.LineTo(bcx+100,bcy);
		dc.MoveTo(900,900); dc.LineTo(1100,900);
		dc.TextOutA(200,70,itoa(N,strtemp,10));
		dc.TextOutA(270,70,itoa(M,strtemp,10));
		dc.TextOutA(340,70,itoa(epsi_count,strtemp,10));
		switch (epsi_umka)
		{
		case 1: dc.TextOutA(410,70,"*2"); break;
		case 2: dc.TextOutA(410,70,"/2"); break;
		case 3: dc.TextOutA(410,70,"/4"); break;
		case 4: dc.TextOutA(410,70,"/10"); break;
		case 5: dc.TextOutA(410,70,"/100"); break;
		case 6: dc.TextOutA(410,70,"=0"); break;
		}
		dc.TextOutA(480,70,itoa(tm_show,strtemp,10));
		dc.TextOutA(550,70,itoa(noise,strtemp,10));
		dc.TextOutA(800,70,itoa(hoffset_fazacf,strtemp,10));
		dc.TextOutA(100,120,itoa(epsi,str2,10));
		//dc.TextOutA(1200,70,itoa(refresh_each_time,strtemp,10));
		dc.TextOutA(1350,70,itoa(scale_rihist,strtemp,10));
		dc.TextOutA(1350,95,itoa(scale_rihist_hor, strtemp, 10));
		flc = 0; //f line count
		csc = 0;
		ifs.clear();
		str2_cfout_make();
		dc.TextOutA(100,100,str2); 
		ifs.open(str2, std::ifstream::in);
		//ifs.open("F:/main/D/! имморталистика/! BRAIN/cerebellum/Projects 5 мая 2014/cerebellum 333112/cerebellum 333/02072014-MNvar2-600/cfout_M4_step0_sin0_noise0_eligc2_epsi3_umka3_tm5018000.txt", std::ifstream::in);
		if (!ifs) 
		{
			if (create_file!=1) 
			{
				dc.SetTextColor(RGB(255,0,0)); 
				dc.TextOutA(100,120,"file doesn't exist. Create?"); 
				dc.SetTextColor(RGB(0,0,0));
			}
			if ((create_file==1)||(refresh_each_time==1)) main_calculation(); 
			create_file=0;
		}
		ifs.close(); ifs.clear();
		str2_cfout_make();
		ifs.open(str2, std::ifstream::in);
		if (ifs)
		{
		ifs >> temp >> temp >> sum >> fazacf;
		for (i=0;i<300;i++) interspike[i]=0;
		while((!ifs.eof())&&(flc<30000))
		{ 			
			//dc.MoveTo(200,200); dc.LineTo(300,300+flc);
			ifs >> temp >> temp >> sum >> fazacf;
			flc++;
			if (flc>offset)
			{
				dc.MoveTo(100+(flc-offset)/scale-1,300-5.*(sum_prev-120.)); dc.LineTo(100+(flc-offset)/scale,300-5.*(sum-120.));
				if (!(flc%1000)) {dc.MoveTo(100+flc/scale,280); dc.LineTo(100+flc/scale,320);};
				if ((fazacf-1!=fazacf_prev) && (fazacf - fazacf_prev != 100-2222))
					{dc.MoveTo(100+(flc-offset)/scale,300); dc.LineTo(100+(flc-offset)/scale,150);};
			}
			if ((fazacf-1!=fazacf_prev) && (fazacf - fazacf_prev != 100-2222))
			{
				//if (csc>hoffset_fazacf) 
				{
						dc.MoveTo(700+5*(csc-hoffset_fazacf),650); 
						dc.LineTo(700+5*(csc-hoffset_fazacf),650-(fazacf_prev-offset_fazacf)*scale_fazacf);
						if ((show_num)&&(!(csc%3))) dc.TextOutA (200+5*(csc-hoffset_fazacf),520+20*(csc%2),itoa(fazacf_prev,strtemp,10));
				}
				csc++;
				dc.SetPixel(bcx,bcy,RGB(255,0,0)); dc.SetPixel(bcx+1,bcy+1,RGB(255,0,0)); dc.SetPixel(bcx,bcy+1,RGB(255,0,0)); dc.SetPixel(bcx+1,bcy,RGB(255,0,0));
				if ((show_corresp) && (!(csc%15))) {dc.MoveTo(200+5*(csc-hoffset_fazacf)+2,500); dc.LineTo(100+(flc-offset)/scale-1,300);};
				interspike_prev = fazacf_prev;
				if (csc<300) interspike[csc] = fazacf_prev; //histogram shows only first 300 intersp.interv. - for normal cases it's enough
			}
			fazacf_prev=fazacf; sum_prev=sum;
		} 
		/*binn=10; for (i=0;i<50;i++) bin[i]=0;
		for (i=0;i<50;i++) for (csc=0;csc<300;csc++) if ((interspike[csc]>i*binn)&&(interspike[csc]<i*binn+binn)) bin[i]++;
		for (i=0;i<50;i++) {dc.MoveTo(900+3*i,900); dc.LineTo(900+3*i,900-bin[i]);};*/
		ifs.close();
		flc = 0; //f line count
		csc = 0;
		ifs.clear();
		str2_timcs_make();
		ifs.open(str2, std::ifstream::in);
		if (!ifs) dc.TextOutA(10,20,"failed open file");
		ifs >> temp >> interspike[0];
		ifs >> temp >> interspike[0];
		for (i=0;i<3000;i++) interspike[i]=0;
		csc=0; interspike_sum_max = 0;
		while((!ifs.eof())&&(flc<30000))
		{ 			
			//dc.MoveTo(200,200); dc.LineTo(300,300+flc);
			ifs >> temp >> interspike[csc];
			flc++;
			interspike_sum_max += interspike[csc];
				dc.SetPixel(bcx+interspike[(csc-1)],bcy-interspike[csc],RGB(0,0,0));//Puankare graph
				dc.SetPixel(bcx,bcy,RGB(255,0,0)); dc.SetPixel(bcx+1,bcy+1,RGB(255,0,0)); dc.SetPixel(bcx,bcy+1,RGB(255,0,0)); dc.SetPixel(bcx+1,bcy,RGB(255,0,0));
		csc++;
		} 
		interspike_sum = 0;
		dc.MoveTo(1200,800); dc.LineTo(1900,800);
		for (i=0;i<csc;i++) //graph to the right of histogram
		{
			i_float = i; csc_float = csc;
			interspike_sum += interspike[i];
			if (i==0) dc.MoveTo(1200+scale_rihist_hor*(i-rihist_offset)/50.,800);
			if (i>rihist_offset)
			{
				//dc.MoveTo(1200+scale_rihist_hor*(i-rihist_offset)/50.,800); 
				dc.LineTo(1200+scale_rihist_hor*(i-rihist_offset)/50.,j=800-300*scale_rihist*(interspike_sum/interspike_sum_max-i_float/(csc_float-1.)));
			}
		}
		binn=10; for (i=0;i<50;i++) bin[i]=0; bin_integral=0;
		for (i=0;i<50;i++) for (csc=0;csc<3000;csc++) if ((interspike[csc]>i*binn)&&(interspike[csc]<i*binn+binn)) {bin[i]++; bin_integral++;};
		for (i=1;i<50;i++) {dc.MoveTo(900+3*i,900); dc.LineTo(900+3*i,900-bin[i]/(bin_integral/2000.));}; //histogram of interspike intervals
		ifs.close();
		flc = 0; //f line count
		csc = 0;
		ifs.clear();
		str2_sgm_make();
		ifs.open(str2, std::ifstream::in);
		if (!ifs) dc.TextOutA(10,20,"failed open file");
		i=0;
		while(i<1000)
		{
			for (int j=0;j<min(N,12000);j++) ifs >> bufsgm[i][j];
			i++;
		}
		ifs.seekg(ios::beg);
		if (show_red) for (int k=0; k<i;k++) 
		{
			lplot_sum=0;
			//dc.TextOutA(200+k,180,str2);
			for (int j=0;j<min(N,12000);j++) 
			{
				lplot_sum += abs(bufsgm[k][j]-bufsgm[0][j]);
			}
			dc.SetPixel(200+k,920,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,921,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,922,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,923,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,924,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,925,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,926,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,927,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,928,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,929,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,930,RGB(lplot_sum*20000.*rgb_scale,0,0));
			dc.SetPixel(200+k,931,RGB(lplot_sum*20000.*rgb_scale,0,0));
		}

		} //from "if (!ifs) "
	//********************************************************end
	CDialog::OnPaint();
	}
}

// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CProbeDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CProbeDlg::OnButton1() //"Go on" button
{
create_file = 1;
Invalidate(TRUE);
CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton2()
{
	switch (N)
	{
		case 10: N=20; M=10; break;
		case 20: N=50; break;
		case 50: N=100; break;
		case 100: N=300; break;
		case 300: N=600; break;
		case 600: N=600; break;
	}
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton3()
{
	switch (N)
	{
		case 10: N=10; break;
		case 20: N=10; break;
		case 50: N=20; break;
		case 100: N=50; break;
		case 300: N=100; break;
		case 600: N=300; break;
	}
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton4()
{
	if (N!=10)
	{
	switch (M)
	{
		case 4: if (N!=600) M=10; if (N==600) M=12; break;
		case 10: M=20; break;
		case 12: M=20; break; //only when N==600
		case 20: M=40; break;
		case 40: M=80; break;
		case 80: M=125; break;
		case 125: M=200; break;
		case 200: M=250; break;
		case 250: M=400; break;
		case 400: M=500; break;
		case 500: M=500; break;
	}
	}
	if (N==10)
	{
	switch (M)
	{
		case 2: M=4; break;
		case 4: M=5; break;
		case 5: M=8; break;
		case 8: M=10; break;
		case 10: M=20; break;
		case 20: M=50; break;
		case 50: M=50; break;
	}
	}
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton5()
{
	if (N!=10) 
	{
	switch (M)
	{
		case 4: M=4; break;
		case 10: M=4; break;
		case 12: M=4; break; //12 - only when N==600
		case 20: if (N!=600) M=10; if (N==600) M=12; break;
		case 40: M=20; break;
		case 80: M=40; break;
		case 125: M=80; break;
		case 200: M=125; break;
		case 250: M=200; break;
		case 400: M=250; break;
		case 500: M=400; break;
	}
	}
	if (N==10)
	{
	switch (M)
	{
	case 2: M=2; break;
	case 4: M=2; break;
	case 5: M=4; break;
	case 8: M=5; break;
	case 10: M=8; break;
	case 20: M=10; break;
	case 50: M=20; break;
	}
	}
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton6()
{
	if (epsi_count<10) epsi_count++;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton7()
{
	if (epsi_count>-5) epsi_count--;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton8()
{
	if (epsi_umka<6) epsi_umka++;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton9()
{
	if (epsi_umka>1) epsi_umka--;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton10()
{
if ((tm_show>6000000)&&(tm_show<Tpred)) tm_show+=2000000;
else
	{
	switch (tm_show)
	{
	case 1998000: tm_show=2018000; break;
	case 2018000: tm_show=2998000; break;
	case 2998000: tm_show=3018000; break;
	case 3018000: tm_show=4998000; break;
	case 4998000: tm_show=5018000; break;
	case 5018000: tm_show=6000000; break;
	case 6000000: if (Tpred==6000000) tm_show=6000000; if (Tpred!=6000000) tm_show=8000000; break;
	}
}
Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton11()
{
if ((tm_show>6000000)&&(tm_show>6000000)) tm_show-=2000000;
else
	{
	switch (tm_show)
	{
	case 1998000: tm_show=1998000; break;
	case 2018000: tm_show=1998000; break;
	case 2998000: tm_show=2018000; break;
	case 3018000: tm_show=2998000; break;
	case 4998000: tm_show=3018000; break;
	case 5018000: tm_show=4998000; break;
	case 6000000: tm_show=5018000; break;
	}
}
Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton12()
{
	if (noise<10) noise++;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton13()
{
	scale = scale*2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}
void CProbeDlg::OnBnClickedButton14()
{
	scale = scale/2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton15()
{
	offset+=(int)(100*scale);
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton16()
{
	if (offset>0) offset-=(int)(100*scale);
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton17()
{
	offset_fazacf-=20/scale_fazacf;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton18()
{
	offset_fazacf+=20/scale_fazacf;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton19()
{
	scale_fazacf*=1.5;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton20()
{
	scale_fazacf/=1.5;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton21()
{
	hoffset_fazacf+=5;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton22()
{
	hoffset_fazacf-=5;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	hoffset_fazacf = csc/100.*my_slider1.GetPos(); //GetPos returns from 0 to 100 always
	*pResult = 0;
}

void CProbeDlg::OnBnClickedButton23()
{
	Invalidate(TRUE);
	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton24()
{
	show_corresp = 1 - show_corresp;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton25()
{
	show_num = 1 - show_num;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton26()
{
	sinmode = 1 - sinmode;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

int random(int Upper)
{
	return (int) (((float)Upper)*rand() / RAND_MAX); //returns 0..Upper-1 (checked)
}

  void PrFqCS() //Printout of CS frequeces
  {
	  int Nnach,inti;
	  Nnach=TecFrCell-29999;
	  if (Nnach<0) 
		  Nnach=0;
	  if ((fall1=fopen("CSfrqu.txt","wt"))==NULL)
		  printf("File for %s CANNOT be OPEND!!\n");
	  for (iwb=Nnach;iwb<(TecFrCell+1);iwb++)
	  {
		  inti=FqCS[iwb];
		  fprintf(fall1,"%f\n",inti); //vydacha chastot
	  }
	  fclose(fall1);/**/
  }


//*********************************************************************
//*********************************************************************


//*********************************************************************
 void generation() // generation of a sequence of inputs
 {
	 int kb,kc,kg,kf,kr,kd,kq,kl,kn,ky,is,iw,r,l,numw,tgen,usp;
	 vecti vspom,vspom1,vspom2,schtch,vspoma,faz;
	 for (ky=0; ky<N; ky++)
		schtch[ky]=0;
	 int c;                    //initialization of random nmbrs
	 c=1117;// 8654;//                 //as above
	 srand(unsigned (c));      //as above
//Initial vector
	 Refr=5; //Refractoriness
	 Texci=3; //Duration of excitation
	 for (kq=0;kq<N;kq++) //for each GC 
	 {
		faz[kq]=2000;
		vspom[kq]=0;
	 }
//Choice of an initial vector
	 ja=0;
	 while (ja<wt) //30 random components of vspom[] make 1, others are 0
		{
		 r=random(N);
		 if (vspom[r]==0)
			{
			 vspom[r]=1;
			 vspoma[ja]=r;
			 ja=ja+1;
			 faz[r]=0;//ADDED NOW, March 2, 1999!!!
			}
		}
// Vector of Excitation
	 kb=0;
	 kg=0;
	 while (kb<10)
	 {
		 if (vspom[kg]==1)
		 {
			 faz[kg]=0;
			 kb++;
		 }
		 kg++;
	 }
	 kc=1;
	 kb=0;
	 while (kb<10)
	 {
		 if (vspom[kg]==1)
		 {
			 faz[kg]=kc;
			 kb++;
		 }
		 kg++;
	 }
	 kc++;
	 kb=0;
	 while (kb<10)
	 {
		 if (vspom[kg]==1)
		 {
			 faz[kg]=kc;
			 kb++;
		 }
		 kg++;
	 }
	 for (kl=0; kl<N; kl++)
	 {
		 vspom1[kl]=0;
		 vspom2[kl]=0;
	 }
	 for (kl=0; kl<N; kl++)
	 {
		 if (faz[kl]==2)
			 vspom1[kl]=1;
		 if (faz[kl]==1)
			 vspom2[kl]=1;
	 }


     for (kr=0;kr<(M-1);kr++)//ostalis': (M-2) & (M-1)
	 {
		 for (kg=0;kg<N;kg++)
			 if (faz[kg]<Texci)
				 vhod[kr][kg]=1;
			 else
				 vhod[kr][kg]=0;
	     for (kl=0; kl<N; kl++)
			 faz[kl]=faz[kl]+1;
	     kb=0;
	     if (kr!=(M-2))
			 while (kb<10)
			 {
				 r=random(N);
		         if (faz[r]>Refr)
				 {
					 faz[r]=0;
					 kb++;
				 }
			 }
	     if (kr>(M-2*Refr))
			 for (kl=0; kl<N; kl++)
				 if ((vspom1[kl]==1)||(vspom2[kl]==1))
					 if (faz[kl]>=Texci)
						 faz[kl]=Texci;
	 }

//M-2 - first two thirds
	 for (kg=0;kg<N;kg++)
		if (vspom1[kg]==1)
		{
			faz[kg]=0;
		}
	 for (kg=0;kg<N;kg++)
		if (faz[kg]<Texci)
			vhod[M-2][kg]=1;
		else
			vhod[M-2][kg]=0;

		for (kl=0; kl<N; kl++)
		faz[kl]=faz[kl]+1;
//Now - for M-1
	 for (kg=0;kg<N;kg++)
		if (vspom2[kg]==1)
		{
			faz[kg]=0;
		}

	 for (kg=0;kg<N;kg++)
		if (faz[kg]<Texci)
			vhod[M-1][kg]=1;
		else
			vhod[M-1][kg]=0;
		
//Now: transform vhod[i][j] into adresses: vhoda[i][j]			
	 for (kl=0; kl<M; kl++)
	 {
		 for (kg=0; kg<N; kg++)
			 if (vhod[kl][kg]==1)
			 {
				 vhoda[kl][kg]=kg;
			 }
	 }
//That is ALL!!!
  }

void Demo()
{
	int imt;
	fazat=0;
	others=0;
	for (imt=0;imt<3000;imt++)
	{
		 if (random(1000)>500) 
		      iph=1;
		 else
			  iph=0;
		 others+=2*iph-1;
		 if ((others>13)||(others<-12))
			  others+=-(2*iph-1);
		fazat++;
		if (fazat>(M-1))
			fazat=0;
		purk=purk*0.9;
		for (ik=0; ik<wt; ik++)
		  {
			 irab1=vhoda[fazat][ik];
             grphase[irab1]=0;
			 perem=sgm[irab1];
			 purk=purk+perem;
             grocc[irab1]++;
		  }
		bufs[s]=purk;//
		buff[s]=Pfunc[fNo][fazat];
		bufe[s]=OtherNC*others+buff[s];
		s=s+1;
		if (s>(Nnbc-1))
		  s=0;
	}
}

void StDemo()
{
	int imt,ikt;
	fazat=0;
	for (imt=0;imt<3000;imt++)
	{
		fazat++;
//		fazat=1;
		if (fazat>(M-1))
			fazat=0;
//////////////////////////////////
        for (ikt=0;ikt<100;ikt++)
//////////////////////////////////
        {
		purk=purk*0.9;
		for (ik=0; ik<wt; ik++)
		  {
			 irab1=vhoda[fazat][ik];
             grphase[irab1]=0;
			 perem=sgm[irab1];
			 purk=purk+perem;
             grocc[irab1]++;
		  }
		}
//////////////////////////////////
		bufs[s]=purk;//
		buff[s]=Pfunc[fNo][fazat];
		s=s+1;
		if (s>(Nnbc-1))
		  s=0;
	}
}


void other_PC_influence_calculation()
{
		if (random(1000)>500) 
		      iph=1;
		 else
			  iph=0;
		 others+=2*iph-1;
		 if ((others>13)||(others<-12))
			  others+=-(2*iph-1);
		rabs=OtherNC*others+alf;
}
void renewing_of_weights_sampling_point (int tm_given)
{
			if (tm==tm_given)
			 {
				 for (i=0;i<N;i++)
				 {
					 oldsgm[i]=sgm[i];
				 }
			     rab=0;
			     for(i=0;i<N;i++)
				     rab+=oldsgm[i]*oldsgm[i];
			     OLDSQ=rab;/**/
			 }
}
void fprintf_and_printf()
{
	if (allowfprintf)
	{
				 make_str2_cfout_fprintf();
				 printf ("%s\n",str2);
				 if ((fall=fopen(str2,"wt"))==NULL)
					 printf("File for %s CANNOT be OPEND!!\n");
				 cons=0;
				 if ((tm==Tepsioff)) printf("\n-211- eligc2prev[%d]=%6.2f",qq,eligc2[2]);
				 for (iaa=0;iaa<Nnbc;iaa++) //Nnbc = buffer size
				 {
					 fprintf(fall,"%6.2f\t",(bufs[iaa])); //vydacha vnutri cikla
					 fprintf(fall,"%6.2f\t%6.2f\t%d",(buff[iaa]),bufs[iaa]+buff[iaa],bufe[iaa]); //vydacha vnutri cikla
					 fprintf(fall,"\n");
				 }
				 fprintf(fall,"%8.6f\t%8.6f\t%8.6f\t%8.6f\n",epsi_basal,epsi_basal,epsi_basal,epsi_basal);
				 if (allowexit1==0) fprintf(fall,"%d\t%d\t%d\t%d",tm_ust1,tm_ust1,tm_ust1,tm_ust1);
				 if (allowexit1==1) fprintf(fall,"%d\t%d\t%d\t%d",tm_ust2,tm_ust2,tm_ust2,tm_ust2);
				 fclose(fall);
				
				make_str2_timcs_fprintf(); 
				 printf ("%s\n",str2);
				if ((timcs_file=fopen(str2,"wt"))==NULL)
				 printf("File for %s CANNOT be OPEND!!\n");
				if ((tm==Tepsioff)) printf("\n-213- eligc2prev[%d]=%6.2f",qq,eligc2[2]);
				 for (timcs_i=0;timcs_i<Ncs;timcs_i++)
				 {
					 fprintf(timcs_file,"%6.0f\t%6.0f\n",float(timcs[timcs_i]),float(timcs[timcs_i]-timcs[(timcs_i-1)]));
				 }
				 fclose(timcs_file);

				 make_str2_sgm_fprintf();
				if ((sgm_file=fopen(str2,"wt"))==NULL)
				 printf("File for %s CANNOT be OPEND!!\n");
				for (int i = 0; i<min(Ncs,20000); i++)
				{
					for (int j=0; j<min(N,1200); j++)
					{
						fprintf(sgm_file,"%9.7f\t",bufsgm[i][j]);
					}
					fprintf(sgm_file,"\n\n");
				}

				 fclose(sgm_file);
					printf("Number of complex spikes = %d, tm=%d, epsi=%9.7f\n",Ncs,tm,epsi);
					printf("sum=%6.2f   %6.2f   %6.2f, CS_interval=%6.2f, purk=%6.2f\n",bufs[25]+buff[25],bufs[50]+buff[50],bufs[75]+buff[75],float(Tmsp*2)/(Ncs-Ncs_prev-0.01),purk); 
					Ncs_prev=Ncs;
		if (progmode==2) 
		{
			if ((tm==Tpred-1000)&&(allowexit1==1)&&(noise==noise_max-1)&&(epsi_umka==umka_max-1)) allowexit2=1;
			else if ((tm==Tpred-1000)&&(noise==noise_max-1)&&(epsi_umka==umka_max-1))
			{
				allowexit1=1; printf("\n\nallowexit1=0!"); Sleep(1000);
				epsi_basal=epsi_basal_first;  
				epsi_count=2;
				allowfprintf=0;
			}
		}
	}
}
void phase_control()
{
		fazat++;//Circular Mode
		 if (fazat>(M-1))//Circulation!!!
			 fazat=0; //fazat-M;           
		 if (fazat<0)//Circulation!!!
			 fazat=M-1;//fazat+M;
}
void unknown_function_1 ()
{
		if ((tm/1000)*1000==tm)//once a thousand times (checked) 
			{
              rab=0;
		      for(i=0;i<N;i++)
				  rab+=sgm[i]*sgm[i];
		      SQ=rab;
              rab=0; 
		      for (i=0;i<N;i++)
				  rab+=(sgm[i]-oldsgm[i])*(sgm[i]-oldsgm[i]);
		      DST=rab/(SQ*OLDSQ);
//              bufe[q]=DST;
              q++;
		      if (q==Nnbc)
				  q=0;
			}
}
void initialization_of_everything()
{
//%%%%%%%%%%%%%%%%%%% CONSISTS OF: 
//%%%%%%%%%%%%%%%%%%% (1) phase function initialization
//%%%%%%%%%%%%%%%%%%% (2) eligibility initialization
//%%%%%%%%%%%%%%%%%%% (3) synaptic plasticity initialization
	Tjump1=2000000;
	Tjump2=3000000;
	Tepsioff=35000000;
	fNo=0; Nnbc = 2000; Nnbc_rab=20000;
	base=60;//80
	fi_cylinder=180.; r_cylinder=noise*noise/500.;
	purk_mult = 0.9; //pow(0.9,(tau/10.));
	intgr_mult = 0.6; //pow(0.6,(tau/10.));
	fazacf=-1;
	if (sinmode==0)
	{
	for (i=0;i<M;i++)  //Pfunc[0][...]  - constant
     {
		 if (i<(M/2))
			 Pfunc[0][i]=base-step;//100;//80;//195
		 else
			 Pfunc[0][i]=base-step;
	 }
	}
	else 
	{
	for (i=0;i<M;i++) Pfunc[0][i]=base+floor(step*sin(2.*3.141592654*i/M));//100;//80;//195
	}
	 for (i=0;i<M;i++)  //Pfunc[1][...]  - saw-tooth
		 if (i<(M/2))
			 Pfunc[1][i]=(i/(M/100.))/2+(67/(M/100.));//120+M/2+i
		 else
			 Pfunc[1][i]=(117/(M/100.))-(i/(M/100.))/2; //120+M-i
	 for (i=0;i<M;i++) printf("%d ",Pfunc[1][i]);
	 printf("\n333333333\n");
	 fNo=0;
//%%%%%%%%%%%%%%%%%%% (2) eligibility initialization
	 for (i=0; i<N; i++)  //Zeroing eligc[*]
		 eligc[i]=eligc1[i]=eligc2[i]=eligc3[i]=0;
//%%%%%%%%%%%%%%%%%%% (2) eligibility initialization: END.
//%%%%%%%%%%%%%%%%%%% (3) synaptic plasticity initialization
	 hhi[0]=-26;
	 hhi[1]=-27;
	 hhi[2]=-27;
	 hhi[3]=-26;
	 hhi[4]=-26;
	// for (int iii=0; iii<(5*10/tau);iii++) hhi[iii]=-26.4;
	 printf("\n6666666666\n");
	for (i=5;i<100000;i++) // for (i=(int)(5*10/tau);i<100000;i++)
		hhi[i]=1;
	printf("\n7777777777\n");
//	 Tequ=(Tmin2-Tmin1)*(47+1);
//%%%%%%%%%%%%%%%%%%% END OF THE FIRST SECTION
//%%%%%%%%%%%%%%%%%%% END OF THE FIRST SECTION
//	 tminus = 30;   // the model's TMINUS
	 FphaseVal=-29; // constant of changes in the first phase
	 SphaseVal=1;   // constant of changes in the second phase
	 generation2(); // generation of a sequence of inputs___________________________________________
	 btt=1;
	 intgr=0;printf("\n555555555555\n");
	 purk=0.;
	 fazat=0;
	 fazacf=0;
	 Ipred=300;//300//50000;
	 Tmsp=5000;//90000;//100000;//1000000;//parameter for control of program 
	 alfama=100;
	 for (ic=0; ic<N; ic++) // weights initialization PC#0//ic - index of granule cell
	 {
         sgm[ic]=10.*random(wt+1)/float(wt); //weight of granule cells synapsw on PC 
	 }
	 for (ic=0; ic<TavP; ic++) //initialization //TavP
	 {
		 ttek[ic]=9000000;//tekushshee
	 }
	 tm=0;
	 printf("\n4444444444444\n");
	 q=0; ///////////////////////////////////////////////
	 s=0;
	 cht=0;
	 Ncs_pred=1;
	 tmg=0; //Global Time
	 Ncs=0; Ncs_prev=0;
	 Nss=0;
	 fazat=0;
	 intpurk=0;
	 ntek=0;
	 Nav=0;
	 Nav1=0;
	 TavPstoh=TavP;
	 ss=0;
	 cs=0;
	 bdip=0;
	 tekcs[0]=400000000;
	 fazacf=1000;
	 fazat=0;
	 razr=1;
	 Tf=1000;
//Below: initialization of oldsgm[i] and OLDSQ
	 for (i=0;i<N;i++)
	 {
		 oldsgm[i]=sgm[i]; //oldsgm - old weight of GC=>PC
	 }
	 rab=0;
	 
	 for(i=0;i<N;i++)
		 rab+=(oldsgm[i])*(oldsgm[i]);
	 OLDSQ=rab; //squares of oldsgm
     others=0; //currently there are no other Nuclear Cells
	 OtherNC=0;//0.3;//1;
     
}
void general_randomization(int rndmz)
{
	 srand(unsigned (rndmz));   //GENERAL RANDOMIZATION
}
void protecting_overflow_of_grphase()
{
		for (qq=0;qq<N;qq++) if (grocc[qq]==0) //It is necessary to protect overflow of grphase
		{
			grphase[qq]=1000;
			sgm[qq]=0;
		}
}
void bufferization_of_purk_and_alf()
{
	bufs[s]=purk;
	buff[s]=alf;
	bufe[s]=fazacf;   
	for (int i=0;i<N;i++) bufsgm[s][i]=sgm[i];
/*	if ((abs(bufe[s]-bufe[(s-1)])>1)&&(fazacf>-1)&&(fazacf!=2222)&&(fazacf!=100)&&(s>0))
	{
		printf("\ntermination of universe in 10 sec\n"); 
		printf ("fazat=%d,fazacf=%d,bufe[s]=%d,tm=%d",fazat,fazacf,bufe[s],tm);
		Sleep(7000);
	};*/
	s++;		  
	if (s>(Nnbc-1)) s=0;
}
void show_gc_activity()
{
	int ifactive;
	for (int i=0; i<M; i++)
	{
		printf ("\n\n");
		ifactive=0;
		for (int j=0; j<N; j++)
		{
			printf ("%d  ",vhoda[i][j]); //i = fazat, j = num of GC
			if (!(vhoda[i][j]==0)) ifactive++;
		}
		printf ("\nnumber of active GCs = %d",ifactive);
	}
Sleep(5000);
	for (int j=0; j<N; j++)
	{
		ifactive=0;
		for (int i=0; i<M; i++)
		{
			//printf ("%d  ",vhoda[i][j]); //i = fazat, j = num of GC
			if (!(vhod[i][j]==0)) ifactive++;
		}
		printf ("\nthis GC has been active  %d times",ifactive);
	}
}
void show_gc_to_pc_weights()
{
	if ((tm%1000000<10)&&(tm%1000000>0)) 
	{
	for (int i=0; i<300; i++) printf ("%6.2f\t%d\t%10.8f\n",sgm[i],vhoda[fazat][i],epsi*sgm[i]*hhi[fazacf]*vhod[fazat][i]);
	printf("\nfazat =%d, fazacf=%d,liana integral=%6.2f of %d\n",fazat,fazacf,intgr,Ipred);
	printf("\ntm=%d, purk=%6.2f",tm,purk);
	Sleep(1000);
	}
}
void generation2()
{
	int i,j,k,fazat,i_rand,exc_left[Nmax],tupic, Nra,Mra; //exc_left[N] - how much excitations left for neurons
	int sum_exc_left, partial_sum;
	fazat=i=j=k=tupic=0;
	for (i=0;i<M;i++) for (j=0;j<N;j++) vhod[i][j]=0;
	for (i=0;i<N;i++) exc_left[i]=(wt*M)/N; //WARNING - may turn to be NOT INTEGER - check it
	if (!((wt*M)%N==0)) {printf("\nWarning or error:(wt*M)\%N is not zero"); Sleep(10000);};
	sum_exc_left = wt*M;
printf("\n888888888\n"); Sleep(1000);
	for (fazat=0;fazat<M;fazat++)
	{
		printf("\nfazat=%d",fazat);
	i=0;
	while ((i<wt)&&(sum_exc_left>1))//"+1" to evade problems of 2/0=>1/1
	{		
		i_rand = random(sum_exc_left); 
		partial_sum = -1;  j= -1;
		while (i_rand > partial_sum)
		{
			j++;
			partial_sum += exc_left[j];
		}
		if (vhod[fazat][j]==0) 
		{
			vhod[fazat][j]=1;
			exc_left[j]--;
			sum_exc_left--;
			i++; tupic=0;
		}
		tupic++;
		if (tupic>10000) //it means that vhod[fazat][j] is already 1 but nobody else can take this increase
		{
			printf("\ntupic>10000");
			Nra=random(N);
			Mra=random(fazat);
			if ((vhod[fazat][Nra]==0)&&(vhod[Mra][Nra]==1)&&(vhod[Mra][j]==0))
			{
				vhod[fazat][Nra]=1;
				vhod[Mra][Nra]=0;
				vhod[Mra][j]=1;
				i++; tupic=0;
				exc_left[j]--;
				sum_exc_left--;
				printf("\ni=%d, fazat=%d (there was tupic)\n",i,fazat);Sleep(500);
			}
		}
		/*if ((fazat==M-1)){ 
		printf("\nsum_exc_left=%d, i_rand=%d,i=%d<wt,fazat=%d,j=%d\n",sum_exc_left,i_rand,i,fazat,j);
		for (int ijk=1;ijk<N;ijk++) printf("%d/%d ",exc_left[ijk],vhod[fazat][ijk]);
		Sleep(1000);};*/
	}
	} //end of "fazat" cycle
	/*while (sum_exc_left>-1)
	{
		i_rand=random(N);
		if (vhod[(M-1)][i_rand]==0)
		{
			vhod[(M-1)][i_rand]=1;
	//		printf("j=%d,fazat=%d  intermittent",j,M-1); Sleep(4000);
			sum_exc_left--;
		}
	};*/
	for (int ijk = 1; ijk<N; ijk++)
	{
		if (exc_left[ijk]>0)
		{
			exc_left[ijk]--;
			vhod[(M-1)][ijk]=1;
		}
	}
printf("\n999999999999\n");
	//Now: transform vhod[i][j] into adresses: vhoda[i][j]			
	 for (i=0; i<M; i++)
	 {
		 int numw=0;
		 for (j=0; j<N; j++)
			 if (vhod[i][j]==1)
			 {
				 vhoda[i][numw]=j;
				 numw++;
			 }
	 }
	//That is ALL!!!

}
void eligibility_as_grphase_initialization()
{
	 for (i=0; i<TW; i++)  //Zeroing eligc[*]
		 eligc[i]=0;
     eligc[1]=1;//Exponential
	 for (i=2; i<1000; i++)        //EXP
	 {                             //EXP
		 rab=0.90;                     //EXP
		 eligc[i]=eligc[i-1]*rab;   //EXP
//		 printf("%f ",eligc[i]);    //EXP
	 }  /**/                         //EXP
     rab=0;//N, normalizing
	 for (i=0; i<TW; i++)//N
		 rab=rab+eligc[i];//N
	 printf("\n Integral of Eligc(tay)=%f\n",rab);//N
	 for (i=0; i<TW; i++)//N
		 eligc[i]=eligc[i]/rab;//N
     rab=0;//N
	 for (i=0; i<TW; i++)//N
		 rab=rab+eligc[i];//N: Normalizing eligc[*]
	 printf("\n Normalized Integral of Eligc(tay)=%f\n",rab);
	 
}

void make_jump1()
{
	for (i=0;i<M;i++)  Pfunc[0][i]=Pfunc[0][i]+20;
}
void make_jump2()
{
	for (i=0;i<M;i++)  Pfunc[0][i]=Pfunc[0][i]-20;
}
void make_epsi_little()
{
	Ncs=1;
	if (epsi_umka==1) epsi=epsi*2.;//2.
	if (epsi_umka==2) epsi=epsi/2.;//3.
	if (epsi_umka==3) epsi=epsi/4.;//7.
	if (epsi_umka==4) epsi=epsi/10.;//1.5
	if (epsi_umka==5) epsi=epsi/100.;//2.
	if (epsi_umka==6) epsi=0.;//1.5
}
void initialize_N_M_wt()
{
	printf("\ninitialNMwt");
	if (Ncount==1)
	{
			Nnb=20000; wt=5;
			N=10;  printf("\nN=10");
			if (Mcount==1) M=2;//{M=4; wt=10;};
			if (Mcount==2) M=4;//{M=10; wt=10;};
			if (Mcount==3) {M=5; wt=4;};//{M=20;};
			if (Mcount==4) M=8;//M=40;
			if (Mcount==5) M=10;//M=80;
			if (Mcount==6) M=20;//{M=125; wt=8;};
			if (Mcount==7) M=50;//M=200;
			//if (Mcount==8) //{M=250;};
			//if (Mcount==9) //M=400;
			if (!(Nnb%M==0)) {printf("Nnb%M != 0, error"); Sleep(100000);};//also, (M*wt)%N has to be 0 [Ctrl+F "WARNING" for proof]
									//also, both wt and M must always be less than N
	}
	if (Ncount==2)
	{
			Nnb=20000; wt=5;
			N=20; 
			if (Mcount==1) {M=4; wt=10;};
			if (Mcount==2) {M=10; wt=10;};
			if (Mcount==3) {M=20;};
			if (Mcount==4) M=40;
			if (Mcount==5) M=80;
			if (Mcount==6) {M=125; wt=8;};
			if (Mcount==7) M=200;
			if (Mcount==8) {M=250;};
			if (Mcount==9) M=400;
			if (!(Nnb%M==0)) {printf("Nnb%M != 0, error"); Sleep(100000);};//also, (M*wt)%N has to be 0 [Ctrl+F "WARNING" for proof]
									//also, both wt and M must always be less than N
	}
	if (Ncount==3)
	{
			Nnb=20000; wt=10;
			N=50; 
			if (Mcount==1) {M=4; wt=25;};
			if (Mcount==2) {M=10; wt=20;};
			if (Mcount==3) {M=20;};
			if (Mcount==4) M=40;
			if (Mcount==5) M=80;
			if (Mcount==6) {M=125; wt=40;};
			if (Mcount==7) M=200;
			if (Mcount==8) {M=250;};
			if (Mcount==9) M=400;
			if (!(Nnb%M==0)) {printf("Nnb%M != 0, error"); Sleep(100000);};//also, (M*wt)%N has to be 0 [Ctrl+F "WARNING" for proof]
									//also, both wt and M must always be less than N
	}
	if (Ncount==4)
	{
			Nnb=20000; wt=10;
			N=100; 
			if (Mcount==1) {M=4; wt=50;};
			if (Mcount==2) {M=10; wt=30;};
			if (Mcount==3) {M=20;};
			if (Mcount==4) M=40;
			if (Mcount==5) M=80;
			if (Mcount==6) {M=125; wt=40;};
			if (Mcount==7) M=200;
			if (Mcount==8) {M=250;};
			if (Mcount==9) M=400;
			if (!(Nnb%M==0)) {printf("Nnb%M != 0, error"); Sleep(100000);};//also, (M*wt)%N has to be 0 [Ctrl+F "WARNING" for proof]
	}
	if (Ncount==5)
	{
			Nnb=20000; wt=30;
			N=300; 
			if (Mcount==1) {M=4; wt=75;};
			if (Mcount==2) {M=10; wt=60;};
			if (Mcount==3) {M=20;};
			if (Mcount==4) M=40;
			if (Mcount==5) M=80;
			if (Mcount==6) {M=125; wt=24;};
			if (Mcount==7) M=200;
			if (Mcount==8) {M=250;};
			if (Mcount==9) M=400;
			if (!(Nnb%M==0)) {printf("Nnb%M != 0, error"); Sleep(100000);};//also, (M*wt)%N has to be 0 [Ctrl+F "WARNING" for proof]
	}
	if (Ncount==6)
	{
			Nnb=20000; wt=30;
			N=600; 
			if (Mcount==1) {M=4; wt=300;};
			if (Mcount==2) {M=12; wt=100;};
			if (Mcount==3) {M=20; wt=120;};
			if (Mcount==4) M=40;
			if (Mcount==5) M=80;
			if (Mcount==6) {M=125; wt=24;};
			if (Mcount==7) M=200;
			if (Mcount==8) {M=250; wt=24;};
			if (Mcount==9) M=400;
			if (!(Nnb%M==0)) {printf("Nnb%M != 0, error"); Sleep(100000);};
	}
}
void initialize_wt()
{
	if (N==10)
	{
			wt=5;
			if (M==5) wt=4;
	}
	if (N==20)
	{
			wt=5;
			if (M==4) wt=10;
			if (M==10) wt=10;
			if (M==125) wt=8;
	}
	if (N==50)
	{
			wt=10;
			if (M==4) wt=25;
			if (M==10) wt=20;
			if (M==125) wt=40;
	}
	if (N==100)
	{
			wt=10;
			if (M==4) wt=50;
			if (M==10) wt=30;
			if (M==125) wt=40;
	}
	if (N==300)
	{
			wt=30;
			if (M==4) wt=75;
			if (M==10) wt=60;
			if (M==125) wt=24;
	}
	if (N==600)
	{
			wt=30;
			if (M==4) wt=300;
			if (M==12) wt=100;
			if (M==20) wt=120;
			if (M==125) wt=24;
			if (M==250) wt=24;
			if (M==500) wt=24;
	}
}
int main_calculation()
	{
		allowfprintf=1;
		
			if (eligc_mode==1) {freq1=0.7; freq2=0.5;}; //below is (exp-exp)^2; it's 0 both when x=0 and x=inf
			if (eligc_mode==2) {freq1=0.9; freq2=0.7;}; //below is (exp-exp)^2; it's 0 both when x=0 and x=inf
			if (eligc_mode==3) {freq1=0.99; freq2=0.7;}; //below is (exp-exp)^2; it's 0 both when x=0 and x=inf
			if (eligc_mode==4) {freq1=0.99; freq2=0.9;}; //below is (exp-exp)^2; it's 0 both when x=0 and x=inf
			if (eligc_mode==0) eligibility_as_grphase_initialization();
			epsi_basal = epsi_basal_first = 0.0000001; //epsi_basal_first is used in other parts of program
//	printf("\ninitialization of everything - begin");
	initialize_wt(); 
	initialization_of_everything();//________________________________________
//	printf("\ninitialization of everything - end");
	Nnb=20000;
//	general_randomization(58654);
//  show_gc_activity(); Sleep(5000);
	printf("\nM=%d",M);
	//epsi = epsi_basal*100000000.;
	epsi = pow(sqrt(10.),epsi_count-1.)*epsi_basal;
	tm_ust1=tm_ust2=0;
	Sleep(500);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 for (tm=0; tm<Tpred; tm++) //main cycle begins 
	 {
		 alf=Pfunc[0][fazat];//  /\/\/\/\/  ("fazat" is const at each iteration of "tm")
		 if (noise>0) 
		 {
			 fi_cylinder+=random(9)-4; if (fi_cylinder>360) fi_cylinder-=360; if (fi_cylinder<0) fi_cylinder+=360;
			 alf=alf+r_cylinder*sin(fi_cylinder*3.141592654/180.);
		 };
// other_PC_influence_calculation(); //currently we have only one purkinje cell		
		 purk=purk*0.9;	
		 for (ik=0; ik<wt; ik++) //updating output of PC
		 {								//wt is number of non-zero components
			  irab1=vhoda[fazat][ik]; //number of active GC if it's active
              grphase[irab1]=0; //so we can to zero it's grphase
			  grocc[irab1]++; //integrating each GC excitations
              purk+=sgm[irab1]; //adding weight of GC=>PC
		 }
		 for (qq=0;qq<N;qq++) //updating weights of GC=>PC 
         {
			 if (eligc_mode==0)
			 sgm[qq]+=float(epsi*hhi[fazacf]*eligc[(grphase[qq])]);//Cerebellar Mode
			 else
			 {
			eligc[qq]=eligc[qq]*0.9+vhod[fazat][qq];//vhod[fazat][qq];//
			eligc1[qq]=eligc1[qq]*(freq1)+vhod[fazat][qq];
			eligc2[qq]=eligc2[qq]*(freq2)+vhod[fazat][qq];
			eligc3[qq]=eligc3[qq]*(freq1+freq2)/2.+vhod[fazat][qq];
			eligc[qq]=eligc1[qq]+eligc2[qq]-2.*eligc3[qq];
			sgm[qq]+=float(epsi*hhi[fazacf]*eligc[qq]);
			 }

//		if (tm%200000<3)	printf("sgm+=%9.8f, fazacf=%d, sgm=%6.2f, eligc=%6.5f,grphase=%d,qq=%d\n",float(epsi*hhi[fazacf]*eligc[(grphase[qq])]), fazacf,sgm[qq],eligc[(grphase[qq])],grphase[qq],qq);
			if (sgm[qq]<0) sgm[qq]=0;
			if (sgm[qq]>5) sgm[qq]=5; //in the beginning, they're initialized up to 10...
			grphase[qq]+=1; if (grphase[qq]>100000) grphase[qq]=1000;
		 }

		if ((tm==Tepsioff)) printf("\n-1- eligc2prev[%d]=%6.2f",qq,eligc2[2]);
//		 show_gc_to_pc_weights(); //weights and their changes
		 intgr=intgr*0.6+alf+purk; //0.6 -standard; //in established regime 300 = 300*0.6 + 120 <=
		 if (intgr>Ipred) //generation of impulse in liana cell
		 {
			  intgr=0.;
			  fazacf=-1;
			  timcs[Ncs]=tm;//Ncs=number of complex spike; timcs=time of CS #Ncs
			  Ncs++;  //if (Ncs==1000000) return 0;
			  if (Ncs>50173) Ncs=1;
		 }
//unknown_function_1(); //it's function of "q" and we don't use "q" now
bufferization_of_purk_and_alf();
protecting_overflow_of_grphase();
		 if (((!(tm%100000)) || (tm>Tepsioff)&&(!(tm%50000))))
		 {
			printf ("tm=%d, tms=%d, Tmsp=%d\n",tm,tms,Tmsp);
			printf("Number of complex spikes = %d, tm=%d, epsi=%9.7f\n",Ncs,tm,epsi);
			printf("sum=%6.2f   %6.2f   %6.2f, purk=%6.2f\n",bufs[25]+buff[25],bufs[50]+buff[50],bufs[75]+buff[75],purk); //*2 for cht
		 }	
		 tms++; if (tms>Tmsp) tms=1;
		 if (tm==Tjump1-22000) {s=tms=1; Tmsp=Nnbc=Nnbc_rab;};
		 if (tm==Tjump1-2000) {fprintf_and_printf(); s=tms=1; Tmsp=Nnbc=Nnbc_rab;};
//		 if (tm==Tjump1) make_jump1();
		 if ((tm>Tjump1)&&(!tm_ust1)) if (abs((purk+alf)-Ipred*0.4)<0.5) 
			{tm_ust1=tm-Tjump1; printf("\n\ntm_ust1=%d\nallowexit1=%d",tm_ust1,allowexit1);Sleep(1000);};
		 if (tm==Tjump1+18000) fprintf_and_printf();
		 if (tm==Tjump2-22000) {s=tms=1; Tmsp=Nnbc=Nnbc_rab;};
		 if (tm==Tjump2-2000) {fprintf_and_printf(); s=tms=1; Tmsp=Nnbc=Nnbc_rab;};
//		 if (tm==Tjump2) make_jump2();
		 if ((tm>Tjump2)&&(!tm_ust2)) if (abs((purk+alf)-Ipred*0.4)<0.5) 
			{tm_ust2=tm-Tjump2; printf("\n\ntm_ust2=%d\nallowexit1=%d",tm_ust2,allowexit1);Sleep(1000);};
		 if (tm==Tjump2+18000) fprintf_and_printf();
		 if (tm==Tepsioff-22000) {s=tms=1; Tmsp=Nnbc=Nnbc_rab;};
		 if (tm==Tepsioff-2000) {fprintf_and_printf(); s=tms=1; Tmsp=Nnbc=Nnbc_rab;};
         if (tm==Tepsioff) {make_epsi_little(); Ipred = Ipred * 1.03;};
		 if (tm==Tepsioff+18000) fprintf_and_printf();
		 if ((!(tm%2000000))) fprintf_and_printf();
		 if ((tm>Tepsioff)&&((tm%2000000)==3)) Ncs=1;
		 if (tm==Tpred) fprintf_and_printf();
		 fazacf++;
		 if (fazacf>2222) fazacf=100; //protecting unhadled violation when fazacf gets 2000000 //not "fazacf=-1" for hhi[fazacf]=1 under any fazacf>5
		 phase_control(); //fazat++ and cycling
	} //tm cycle ends;

	return(0);
}  //end of main

void CProbeDlg::OnBnClickedButton27()
{
	step=step+5;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton28()
{
	if (step>20) step-=5;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}
void CProbeDlg::OnBnClickedButton29()
{
	refresh_each_time = 1 - refresh_each_time;
	Invalidate(TRUE);
	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton30()
{
	if (eligc_mode<5) eligc_mode++;
	if (eligc_mode==5) {freq1=0.7; freq2=0.5;};
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton31()
{
	if (eligc_mode>0) eligc_mode--;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton32()
{
	if (noise>0) noise--;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void str2_cfout_make()
{
			strcpy(str1,"gg");
				 itoa(N,str1,10);
				 strcpy(str2,"I:/viz3/cfout_N");
				 strcat(str2,str1);
				 itoa(M,str1,10);
				 strcat(str2,"_M");
				 strcat(str2,str1);
				 itoa(step,str1,10);
				 strcat(str2,"_step");
				 strcat(str2,str1);
				 itoa(sinmode,str1,10);
				 strcat(str2,"_sin");
				 strcat(str2,str1);
				 itoa(noise,str1,10);
				 strcat(str2,"_noise");
				 strcat(str2,str1);
				 if (eligc_mode!=5)
				 {
				 itoa(eligc_mode,str1,10);
				 strcat(str2,"_eligc");
				 strcat(str2,str1);
				 }
				 else
				 {
				 itoa(100.*freq1,str1,10);
				 strcat(str2,"_1f");
				 strcat(str2,str1);
				 itoa(100.*freq2,str1,10);
				 strcat(str2,"_2f");
				 strcat(str2,str1);
				 }
				 itoa(epsi_count,str1,10);
				 strcat(str2,"_epsi");
				 strcat(str2,str1);
				 itoa(epsi_umka,str1,10);
				 strcat(str2,"_umka");
				 strcat(str2,str1);
				 if (Tpred!=6000000)
				 {
				 itoa(Tpred,str1,10);
				 strcat(str2,"_tpr");
				 strcat(str2,str1);
				 }
				 itoa(tm_show,str1,10);
				 strcat(str2,"_tm");
				 strcat(str2,str1);
				 strcat(str2,".txt");
}
void str2_timcs_make()
{
	 strcpy(str1,"gg");
				 itoa(N,str1,10);
				 strcpy(str2,"I:/viz3/timcs_N");
				 strcat(str2,str1);
				 itoa(M,str1,10);
				 strcat(str2,"_M");
				 strcat(str2,str1);
				 itoa(step,str1,10);
				 strcat(str2,"_step");
				 strcat(str2,str1);
				 itoa(sinmode,str1,10);
				 strcat(str2,"_sin");
				 strcat(str2,str1);
				 itoa(noise,str1,10);
				 strcat(str2,"_noise");
				 strcat(str2,str1);
				 if (eligc_mode!=5)
				 {
				 itoa(eligc_mode,str1,10);
				 strcat(str2,"_eligc");
				 strcat(str2,str1);
				 }
				 else
				 {
				 itoa(100.*freq1,str1,10);
				 strcat(str2,"_1f");
				 strcat(str2,str1);
				 itoa(100.*freq2,str1,10);
				 strcat(str2,"_2f");
				 strcat(str2,str1);
				 }
				 itoa(epsi_count,str1,10);
				 strcat(str2,"_epsi");
				 strcat(str2,str1);
				 itoa(epsi_umka,str1,10);
				 strcat(str2,"_umka");
				 strcat(str2,str1);
				 if (Tpred!=6000000)
				 {
				 itoa(Tpred,str1,10);
				 strcat(str2,"_tpr");
				 strcat(str2,str1);
				 }
				 itoa(tm_show,str1,10);
				 strcat(str2,"_tm");
				 strcat(str2,str1);
				 strcat(str2,".txt");
}
void str2_sgm_make()
{
	strcpy(str1,"gg");
	itoa(N,str1,10);
	strcpy(str2,"I:/viz3/sgm_N");
	strcat(str2,str1);
	itoa(M,str1,10);
	strcat(str2,"_M");
	strcat(str2,str1);
	itoa(step,str1,10);
	strcat(str2,"_step");
	strcat(str2,str1);
	itoa(sinmode,str1,10);
	strcat(str2,"_sin");
	strcat(str2,str1);
	itoa(noise,str1,10);
	strcat(str2,"_noise");
	strcat(str2,str1);
	if (eligc_mode!=5)
				 {
				 itoa(eligc_mode,str1,10);
				 strcat(str2,"_eligc");
				 strcat(str2,str1);
				 }
				 else
				 {
				 itoa(100.*freq1,str1,10);
				 strcat(str2,"_1f");
				 strcat(str2,str1);
				 itoa(100.*freq2,str1,10);
				 strcat(str2,"_2f");
				 strcat(str2,str1);
				 }
	itoa(epsi_count,str1,10);
	strcat(str2,"_epsi");
	strcat(str2,str1);
	itoa(epsi_umka,str1,10);
	strcat(str2,"_umka");
	strcat(str2,str1);
	if (Tpred!=6000000)
	{
		itoa(Tpred,str1,10);
		strcat(str2,"_tpr");
		strcat(str2,str1);
	}
	itoa(tm_show,str1,10);
	strcat(str2,"_tm");
	strcat(str2,str1);
	strcat(str2,".txt");
}
void make_str2_cfout_fprintf()
{
			strcpy(str1,"gg");
				 itoa(N,str1,10);
				 strcpy(str2,"I:/viz3/cfout_N");
				 strcat(str2,str1);
				 itoa(M,str1,10);
				 strcat(str2,"_M");
				 strcat(str2,str1);
				 itoa(step,str1,10);
				 strcat(str2,"_step");
				 strcat(str2,str1);
				 itoa(sinmode,str1,10);
				 strcat(str2,"_sin");
				 strcat(str2,str1);
				 itoa(noise,str1,10);
				 strcat(str2,"_noise");
				 strcat(str2,str1);
				 if (eligc_mode!=5)
				 {
				 itoa(eligc_mode,str1,10);
				 strcat(str2,"_eligc");
				 strcat(str2,str1);
				 }
				 else
				 {
				 itoa(100.*freq1,str1,10);
				 strcat(str2,"_1f");
				 strcat(str2,str1);
				 itoa(100.*freq2,str1,10);
				 strcat(str2,"_2f");
				 strcat(str2,str1);
				 }
				 itoa(epsi_count,str1,10);
				 strcat(str2,"_epsi");
				 strcat(str2,str1);
				 itoa(epsi_umka,str1,10);
				 strcat(str2,"_umka");
				 strcat(str2,str1);
				 if (Tpred!=6000000)
				 {
				 itoa(Tpred,str1,10);
				 strcat(str2,"_tpr");
				 strcat(str2,str1);
				 }
				 itoa(tm,str1,10);
				 strcat(str2,"_tm");
				 strcat(str2,str1);
				 strcat(str2,".txt");
}
void make_str2_timcs_fprintf()
{
	strcpy(str1,"gg");
	itoa(N,str1,10);
	strcpy(str2,"I:/viz3/timcs_N");
	strcat(str2,str1);
	itoa(M,str1,10);
	strcat(str2,"_M");
	strcat(str2,str1);
	itoa(step,str1,10);
	strcat(str2,"_step");
	strcat(str2,str1);
	itoa(sinmode,str1,10);
	strcat(str2,"_sin");
	strcat(str2,str1);
	itoa(noise,str1,10);
	strcat(str2,"_noise");
	strcat(str2,str1);
	if (eligc_mode!=5)
				 {
				 itoa(eligc_mode,str1,10);
				 strcat(str2,"_eligc");
				 strcat(str2,str1);
				 }
				 else
				 {
				 itoa(100.*freq1,str1,10);
				 strcat(str2,"_1f");
				 strcat(str2,str1);
				 itoa(100.*freq2,str1,10);
				 strcat(str2,"_2f");
				 strcat(str2,str1);
				 }
	itoa(epsi_count,str1,10);
	strcat(str2,"_epsi");
	strcat(str2,str1);
	itoa(epsi_umka,str1,10);
	strcat(str2,"_umka");
	strcat(str2,str1);
	if (Tpred!=6000000)
	{
		itoa(Tpred,str1,10);
		strcat(str2,"_tpr");
		strcat(str2,str1);
	}
	itoa(tm,str1,10);
	strcat(str2,"_tm");
	strcat(str2,str1);
	strcat(str2,".txt");
}
void make_str2_sgm_fprintf()
{
	strcpy(str1,"gg");
	itoa(N,str1,10);
	strcpy(str2,"I:/viz3/sgm_N");
	strcat(str2,str1);
	itoa(M,str1,10);
	strcat(str2,"_M");
	strcat(str2,str1);
	itoa(step,str1,10);
	strcat(str2,"_step");
	strcat(str2,str1);
	itoa(sinmode,str1,10);
	strcat(str2,"_sin");
	strcat(str2,str1);
	itoa(noise,str1,10);
	strcat(str2,"_noise");
	strcat(str2,str1);
	if (eligc_mode!=5)
				 {
				 itoa(eligc_mode,str1,10);
				 strcat(str2,"_eligc");
				 strcat(str2,str1);
				 }
				 else
				 {
				 itoa(100.*freq1,str1,10);
				 strcat(str2,"_1f");
				 strcat(str2,str1);
				 itoa(100.*freq2,str1,10);
				 strcat(str2,"_2f");
				 strcat(str2,str1);
				 }
	itoa(epsi_count,str1,10);
	strcat(str2,"_epsi");
	strcat(str2,str1);
	itoa(epsi_umka,str1,10);
	strcat(str2,"_umka");
	strcat(str2,str1);
	if (Tpred!=6000000)
	{
		itoa(Tpred,str1,10);
		strcat(str2,"_tpr");
		strcat(str2,str1);
	}
	itoa(tm,str1,10);
	strcat(str2,"_tm");
	strcat(str2,str1);
	strcat(str2,".txt");
}
void CProbeDlg::OnBnClickedButton33()
{
	Tpred+=1000000;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton34()
{
	if (Tpred>6000000) Tpred-=1000000;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}


void CProbeDlg::OnBnClickedButton35()
{
	freq1+=0.01;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton36()
{
	freq1-=0.01;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton37()
{
	freq2+=0.01;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton38()
{
	freq2-=0.01;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton39()
{
	str2_cfout_make();
	remove(str2);
	Invalidate(TRUE); CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton40()
{
	scale_rihist *=2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton41()
{
	scale_rihist /=2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton42()
{
	show_red = 1 - show_red;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton43()
{
	scale_rihist_hor *=2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton44()
{
	scale_rihist_hor /=2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton45()
{
	rihist_offset += 2000/scale_rihist_hor;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton46()
{
	rihist_offset -= 2000/scale_rihist_hor;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton47()
{
	rgb_scale*=2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}

void CProbeDlg::OnBnClickedButton48()
{
	rgb_scale/=2.;
	Invalidate(TRUE);	CProbeDlg::OnPaint();
}
