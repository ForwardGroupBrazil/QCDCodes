#ifndef QCDEventHdr_h
#define QCDEventHdr_h

class QCDEventHdr 
{
    public:
      //------------ Constructor ------------------------------
      QCDEventHdr() { mRun = 0;}
      //------------ Destructor -------------------------------
      ~QCDEventHdr() {}
      void setRun(int fRun)                                  {mRun   = fRun;}
      void setEvt(int fEvt)                                  {mEvent = fEvt;}
      void setLumi(int fLumi)                                {mLumi  = fLumi;}
      void setBunch(int fBunch)                              {mBunch = fBunch;}
      void setPthat(float fPthat)                            {mPthat = fPthat;}
      void setWeight(float fWeight)                          {mWeight = fWeight;} 
      void setRho(float fCaloRho, float fPFRho)              {mCaloRho = fCaloRho; mPFRho = fPFRho;}
      void setVertices(int fNVtx, int fNVtxGood)             {mNVtx = fNVtx; mNVtxGood = fNVtxGood;}
      void setPV(bool fIsPVgood, float fndof, float fx, float fy, float fz) {mIsPVgood = fIsPVgood; mPVndof = fndof; mPVx = fx; mPVy = fy; mPVz = fz;}
      void setHCALNoise(bool fLoose, bool fTight)            {mLooseHCALNoise = fLoose; mTightHCALNoise = fTight;}
      void setPU(int fNBX, int fOOTPUEarly, int fOOTPULate, int fINTPU) {mNBX = fNBX; mOOTPUEarly = fOOTPUEarly; mOOTPULate = fOOTPULate; mINTPU = fINTPU;}
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
      bool isPVgood()       const {return mIsPVgood;}
      bool looseHCALNoise() const {return mLooseHCALNoise;}
      bool tightHCALNoise() const {return mTightHCALNoise;}
      float PVndof()        const {return mPVndof;} 
      float PVx()           const {return mPVx;}
      float PVy()           const {return mPVy;}
      float PVz()           const {return mPVz;}
      float pthat()         const {return mPthat;}
      float weight()        const {return mWeight;} 
      float caloRho()       const {return mCaloRho;} 
      float pfRho()         const {return mPFRho;} 
      private:
        bool mIsPVgood; 
        bool mLooseHCALNoise;
        bool mTightHCALNoise;
        int mRun;
        int mEvent; 
        int mLumi;
        int mBunch;
        int mNVtx;
        int mNVtxGood;
        int mOOTPUEarly;
        int mOOTPULate;
        int mINTPU;
        int mNBX; 
        float mPVndof;
        float mPVx;
        float mPVy;
        float mPVz;
        float mPthat;
        float mWeight;
        float mCaloRho;
        float mPFRho;
};
#endif
