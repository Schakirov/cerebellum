// probe.h : main header file for the PROBE application
//

#if !defined(AFX_PROBE_H__EBFB7F3C_B140_4589_81F2_284C7C6DBA47__INCLUDED_)
#define AFX_PROBE_H__EBFB7F3C_B140_4589_81F2_284C7C6DBA47__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CProbeApp:
// See probe.cpp for the implementation of this class
//

class CProbeApp : public CWinApp
{
public:
	CProbeApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CProbeApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CProbeApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PROBE_H__EBFB7F3C_B140_4589_81F2_284C7C6DBA47__INCLUDED_)
