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
      void setPrescales(int fPreL1, int fPreHLT)             {mPreL1 = fPreL1; mPreHLT = fPreHLT;}
      void setVertices(int fNVtx, int fNVtxGood)             {mNVtx = fNVtx; mNVtxGood = fNVtxGood;}
      void setPV(bool fIsPVgood, float fndof, float fx, float fy, float fz) {mIsPVgood = fIsPVgood; mPVndof = fndof; mPVx = fx; mPVy = fy; mPVz = fz;}
      int runNo()     {return mRun;} 
      int event()     {return mEvent;} 
      int lumi()      {return mLumi;}
      int bunch()     {return mBunch;}
      int preL1()     {return mPreL1;}
      int preHLT()    {return mPreHLT;}
      int nVtx()      {return mNVtx;}
      int nVtxGood()  {return mNVtxGood;}
      bool isPVgood() {return mIsPVgood;}
      float PVndof()  {return mPVndof;} 
      float PVx()     {return mPVx;}
      float PVy()     {return mPVy;}
      float PVz()     {return mPVz;}
      float pthat()   {return mPthat;}
      float weight()  {return mWeight;} 
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
        float mPVndof;
        float mPVx;
        float mPVy;
        float mPVz;
        float mPthat;
        float mWeight;
};
#endif
