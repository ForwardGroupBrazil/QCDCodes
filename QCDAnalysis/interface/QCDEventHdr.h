#ifndef QCDEventHdr_h
#define QCDEventHdr_h

class QCDEventHdr 
{
    public:
      //------------ Constructor ------------------------------
      QCDEventHdr() { mRun = 0;}
      //------------ Destructor -------------------------------
      ~QCDEventHdr() {}
      //------------ Set methods ------------------------------
      void setRun(int fRun)                                  {mRun   = fRun;}
      void setEvt(int fEvt)                                  {mEvent = fEvt;}
      void setLumi(int fLumi)                                {mLumi  = fLumi;}
      void setBunch(int fBunch)                              {mBunch = fBunch;}
      void setPthat(float fPthat)                            {mPthat = fPthat;}
      void setWeight(float fWeight)                          {mWeight = fWeight;} 
      void setRho(float fCaloRho, float fPFRho)              {mCaloRho = fCaloRho; mPFRho = fPFRho;}
      void setVertices(int fNVtx, int fNVtxGood)             {mNVtx = fNVtx; mNVtxGood = fNVtxGood;}
      void setPV(bool fIsPVgood, float fndof, float fx, float fy, float fz) {mIsPVgood = fIsPVgood; mPVndof = fndof; mPVx = fx; mPVy = fy; mPVz = fz;}
      void setBS(float fBSx, float fBSy, float fBSz) {mBSx = fBSx; mBSy = fBSy; mBSz = fBSz;}
      void setHCALNoise(bool fNoise) {mHCALNoise = fNoise;}
      void setPU(int fNBX, int fOOTPUEarly, int fOOTPULate, int fINTPU) {mNBX = fNBX; mOOTPUEarly = fOOTPUEarly; mOOTPULate = fOOTPULate; mINTPU = fINTPU;}
      void setTrPu(float fTrPu) {mTrPu = fTrPu;} // setting the true PU 
      void setXsec(double fXsec) {mmXsec = fXsec;};
      void setInitialPartons(float finitial1, float finitial2) {mInitial_parton1 = finitial1; mInitial_parton2 = finitial2;}
      //------------ Get methods ------------------------------
      int runNo()           const {return mRun;} 
      int event()           const {return mEvent;} 
      int lumi()            const {return mLumi;}
      int bunch()           const {return mBunch;}
      int nVtx()            const {return mNVtx;}
      int nVtxGood()        const {return mNVtxGood;}
      int ootpuEarly()      const {return mOOTPUEarly;}
      int ootpuLate()       const {return mOOTPULate;}
      int intpu()           const {return mINTPU;}
      int nbx()             const {return mNBX;} 
      int pu()              const {return mOOTPUEarly+mOOTPULate+mINTPU;}
      float trpu()          const {return mTrPu;} // get method for True number of interaction
      bool isPVgood()       const {return mIsPVgood;}
      bool hcalNoise()      const {return mHCALNoise;}
      float PVndof()        const {return mPVndof;} 
      float PVx()           const {return mPVx;}
      float PVy()           const {return mPVy;}
      float PVz()           const {return mPVz;}
      float BSx()           const {return mBSx;}
      float BSy()           const {return mBSy;}
      float BSz()           const {return mBSz;}
      float pthat()         const {return mPthat;}
      float weight()        const {return mWeight;}
      float caloRho()       const {return mCaloRho;} 
      float pfRho()         const {return mPFRho;} 
      float Xsec()          const {return mmXsec;} 
      float init_parton1()  const {return mInitial_parton1;}
      float init_parton2()  const {return mInitial_parton2;}
      private:
        //---- flag about the PV quality -------------- 
        bool mIsPVgood; 
        //---- flags about the HCAL noise -------------
        bool mHCALNoise;
        //---- run number ----------------------------- 
        int mRun;
        //---- event number ---------------------------
        int mEvent; 
        //---- lumi section ---------------------------
        int mLumi;
        //---- bunch crossing -------------------------
        int mBunch;
        //---- number of reconstructed vertices -------
        int mNVtx;
        //---- number of good reco vertice ------------ 
        int mNVtxGood;
        //---- simulated (early) out-of-time pu -------
        int mOOTPUEarly;
        //---- simulated (late) out-of-time pu --------
        int mOOTPULate;
        //---- simulated in-time pu -------------------
        int mINTPU;
        //---- number of simulated bunch crossings ----
        int mNBX; 
        //---- ndof for the primary vertex ------------
        float mPVndof;
        //---- event weight for pileup reweighting --------//
        float mTrPu;
        //---- position of the primary vertex --------- 
        float mPVx;
        float mPVy;
        float mPVz;
        //---- position of the beam spot --------------
        float mBSx;
        float mBSy;
        float mBSz;
        //---- simulated pthat variable ---------------
        float mPthat;
        //---- simulation weight ----------------------
        float mWeight;
        //---- median CALO pt density -----------------
        float mCaloRho;
        //---- median PF pt density -------------------
        float mPFRho;
        // ---- Xsec
        double mmXsec;
        // --- Initial Partons ------------------------
        float mInitial_parton1;
        float mInitial_parton2;
};
#endif
