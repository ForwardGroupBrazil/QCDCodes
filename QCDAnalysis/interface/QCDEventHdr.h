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
      void setPU(float fPU)                                  {mPU = fPU;}
      void setPrescales(int fPreL1, int fPreHLT)             {mPreL1 = fPreL1; mPreHLT = fPreHLT;}
      void setVertices(int fNVtx, int fNVtxGood)             {mNVtx = fNVtx; mNVtxGood = fNVtxGood;}
      void setPV(bool fIsPVgood, float fndof, float fx, float fy, float fz) {mIsPVgood = fIsPVgood; mPVndof = fndof; mPVx = fx; mPVy = fy; mPVz = fz;}
      int runNo()     const {return mRun;} 
      int event()     const {return mEvent;} 
      int lumi()      const {return mLumi;}
      int bunch()     const {return mBunch;}
      int preL1()     const {return mPreL1;}
      int preHLT()    const {return mPreHLT;}
      int nVtx()      const {return mNVtx;}
      int nVtxGood()  const {return mNVtxGood;}
      int PU()        const {return mPU;}         
      bool isPVgood() const {return mIsPVgood;}
      float PVndof()  const {return mPVndof;} 
      float PVx()     const {return mPVx;}
      float PVy()     const {return mPVy;}
      float PVz()     const {return mPVz;}
      float pthat()   const {return mPthat;}
      float weight()  const {return mWeight;} 
      private:
        bool mIsPVgood; 
        int mRun;
        int mEvent; 
        int mLumi;
        int mBunch;
        int mPreL1;
        int mPreHLT;
        int mNVtx;
        int mNVtxGood;
        int mPU;
        float mPVndof;
        float mPVx;
        float mPVy;
        float mPVz;
        float mPthat;
        float mWeight;
};
#endif
