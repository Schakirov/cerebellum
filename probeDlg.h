// probeDlg.h : header file
//

#include "afxcmn.h"
#if !defined(AFX_PROBEDLG_H__A5AC1271_3A8C_4737_B994_59D6AB8ED86A__INCLUDED_)
#define AFX_PROBEDLG_H__A5AC1271_3A8C_4737_B994_59D6AB8ED86A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CProbeDlg dialog

class CProbeDlg : public CDialog
{
// Construction
public:
	CProbeDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CProbeDlg)
	enum { IDD = IDD_PROBE_DIALOG };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CProbeDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CProbeDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnButton1();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedButton4();
	afx_msg void OnBnClickedButton5();
	afx_msg void OnBnClickedButton6();
	afx_msg void OnBnClickedButton7();
	afx_msg void OnBnClickedButton8();
	afx_msg void OnBnClickedButton9();
	afx_msg void OnBnClickedButton10();
	afx_msg void OnBnClickedButton11();
	afx_msg void OnBnClickedButton12();
	afx_msg void OnBnClickedButton13();
	afx_msg void OnBnClickedButton14();
	afx_msg void OnBnClickedButton15();
	afx_msg void OnBnClickedButton16();
	afx_msg void OnBnClickedButton17();
	afx_msg void OnBnClickedButton18();
	afx_msg void OnBnClickedButton19();
	afx_msg void OnBnClickedButton20();
	afx_msg void OnBnClickedButton21();
	afx_msg void OnBnClickedButton22();
	afx_msg void OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult);
	CSliderCtrl my_slider1;
	afx_msg void OnBnClickedButton23();
	afx_msg void OnBnClickedButton24();
	afx_msg void OnBnClickedButton25();
	afx_msg void OnBnClickedButton26();
	afx_msg void OnBnClickedButton27();
	afx_msg void OnBnClickedButton28();
	afx_msg void OnBnClickedButton29();
	afx_msg void OnBnClickedButton30();
	afx_msg void OnBnClickedButton31();
	afx_msg void OnBnClickedButton32();
	afx_msg void OnBnClickedButton33();
	afx_msg void OnBnClickedButton34();
	afx_msg void OnBnClickedButton35();
	afx_msg void OnBnClickedButton36();
	afx_msg void OnBnClickedButton37();
	afx_msg void OnBnClickedButton38();
	afx_msg void OnBnClickedButton39();
	afx_msg void OnBnClickedButton40();
	afx_msg void OnBnClickedButton41();
	afx_msg void OnBnClickedButton42();
	afx_msg void OnBnClickedButton43();
	afx_msg void OnBnClickedButton44();
	afx_msg void OnBnClickedButton45();
	afx_msg void OnBnClickedButton46();
	afx_msg void OnBnClickedButton47();
	afx_msg void OnBnClickedButton48();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PROBEDLG_H__A5AC1271_3A8C_4737_B994_59D6AB8ED86A__INCLUDED_)
