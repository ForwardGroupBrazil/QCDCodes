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
      void setVertices(int fNVtx, int fNVtxGood)             {mNVtx = fNVtx; mNVtxGood = fNVtxGood;}
      void setPV(bool fIsPVgood, float fndof, float fx, float fy, float fz) {mIsPVgood = fIsPVgood; mPVndof = fndof; mPVx = fx; mPVy = fy; mPVz = fz;}
      void setHCALNoise(bool fLoose, bool fTight)            {mLooseHCALNoise = fLoose; mTightHCALNoise = fTight;}
      int runNo()           const {return mRun;} 
      int event()           const {return mEvent;} 
      int lumi()            const {return mLumi;}
      int bunch()           const {return mBunch;}
      int nVtx()            const {return mNVtx;}
      int nVtxGood()        const {return mNVtxGood;}
      bool isPVgood()       const {return mIsPVgood;}
      bool looseHCALNoise() const {return mLooseHCALNoise;}
      bool tightHCALNoise() const {return mTightHCALNoise;}
      float PVndof()        const {return mPVndof;} 
      float PVx()           const {return mPVx;}
      float PVy()           const {return mPVy;}
      float PVz()           const {return mPVz;}
      float pthat()         const {return mPthat;}
      float weight()        const {return mWeight;} 
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
        float mPVndof;
        float mPVx;
        float mPVy;
        float mPVz;
        float mPthat;
        float mWeight;
};
#endif
