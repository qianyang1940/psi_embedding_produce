////////////////////////////////////////////////////////////////////////////////////////////////////
/*!\fn doEmcEmbedEvent
\author Alexandre Suaide
*/
class StChain;
StChain *chain=0;
//void doEmcEmbedEvent(int nevents = 10,char* file="../st_zerobias_7103008_raw_1110001.event.root",Bool_t print = kTRUE);
//void doEmcEmbedEvent(int nevents = 10,char* file="st_zerobias_7072058_raw_1110001.event.root",Bool_t print = kTRUE)

void doEmcEmbedEvent(int nevents =1000,char* file=
//"./st_physics_adc_12094016_raw_1500004.geant.root"
//"/star/data19/embedding/pp500_production_2011/Psi2sEE_100_20142601/P11id.SL11d_embed/2011/038/12038080/st_physics_adc_12038080_raw_2510001.event.root" 
"/star/data19/embedding/pp500_production_2011/Psi2sEE_100_20142601/P11id.SL11d_embed/2011/081/12081022/st_physics_adc_12081022_raw_3510005.event.root"
, char* outDir="./",Bool_t print = kFALSE)
{

    if(gClassTable->GetID("TTable") < 0) {
      gSystem->Load("libStar");
      gSystem->Load("libPhysics");
    }
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    //gROOT->Macro("LoadLogger.C");
    gROOT->Macro("loadMuDst.C");
    gSystem->Load("StarMagField");
    gSystem->Load("StMagF");
    gSystem->Load("StTpcDb");
    gSystem->Load("StDetectorDbMaker");
    gSystem->Load("StDbUtilities");
    gSystem->Load("StEEmcUtil");
    gSystem->Load("StEmcUtil");
    gSystem->Load("StEEmcDbMaker");
    gSystem->Load("StMcEvent"); 
    gSystem->Load("StMcEventMaker"); 
    gSystem->Load("StDaqLib");
    gSystem->Load("St_g2t.so");
    gSystem->Load("St_geant_Maker.so");
    gSystem->Load("StAssociationMaker");
    gSystem->Load("StMcAnalysisMaker");
    gSystem->Load("StDbBroker");
    gSystem->Load("St_db_Maker");
    gSystem->Load("libgeometry_Tables");
    gSystem->Load("StEmcRawMaker");
    gSystem->Load("StEmcADCtoEMaker");
    gSystem->Load("StPreEclMaker");
    gSystem->Load("StEpcMaker");
    gSystem->Load("StEmcSimulatorMaker");     
    //gSystem->Load("StEmcMixerMaker");
    //gSystem->Load("StEmcTriggerMaker");
    gSystem->Load("StTriggerUtilities");

    gSystem->Load("StEvent");
    gSystem->Load("StEventMaker");

    //my analysis maker
    //gSystem->Load("StCheckMaker");
    gSystem->Load("StMyElectronMaker");
    cout<<"Makers are loadded"<<endl;
     
    // create chain    
    chain = new StChain("bfc");  
    if(print) chain->SetDebug(1);
    else chain->SetDebug(0);

    StIOMaker* io = new StIOMaker("IO");
    io->SetFile(file);
    io->SetIOMode("r"); 
    io->SetBranch("*",0,"0");           //deactivate all branches
    io->SetBranch("eventBranch",0,"r");
    io->SetBranch("geantBranch",0,"r");

    St_db_Maker *db1 = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","$PWD/StarDb");
    //new StEEmcDbMaker("eemcDb");
    //cout<<"DB are loadded"<<endl;
    //db1->SetFlavor("sim", "bprsCalib");

 // StEventMaker* eventReader = new StEventMaker("events","title");
 // eventReader->doPrintMemoryInfo = kFALSE;
 // eventReader->doPrintEventInfo = kFALSE;
  //StMcEventMaker*     mcEventReader = new StMcEventMaker; // Make an instance..
  //mcEventReader->doPrintEventInfo = kFALSE;
  //mcEventReader->doPrintMemoryInfo = kFALSE;
  //mcEventReader->doUseEemc = kTRUE;


    //StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
    //Adam and Matt's solution
    //adc->setEmbeddingMode(kTRUE);//no longer needed
    // this line is important to propagate all the hits into StEvent
    // so, even the pedestals are propagated. In this case
    // the second AdcToEMaker will be responsible for making the
    // cuts (after the simulated hits are embedded)
    //adc->saveAllStEvent(kTRUE);
    //if(!print) adc->setPrint(kFALSE);

    //StEmcPreMixerMaker *preMixer = new StEmcPreMixerMaker("preEmbed");

    StMcEventMaker *mcEvent = new StMcEventMaker();

    StEmcSimulatorMaker *emcSim = new StEmcSimulatorMaker();

    //StEmcMixerMaker *emb = new StEmcMixerMaker();
    // include the next line if you want to embedd all simuated hits
    // even the ones that do not have a hit in the real data
    //emb->setEmbedAll(kTRUE);
    //if(!print) emb->setPrint(kFALSE);

    StEmcADCtoEMaker *adc1 = new StEmcADCtoEMaker("EReadEmbed");      
    adc1->setEmbeddingMode(kTRUE);
    if(!print) adc1->setPrint(kFALSE);

    StPreEclMaker *pre = new StPreEclMaker();
    if(!print) pre->setPrint(kFALSE);

    StEpcMaker *epc = new StEpcMaker();
    if(!print) epc->setPrint(kFALSE);

    StAssociationMaker    *association = new StAssociationMaker();       // TPC association maker
    association->useInTracker();
    //association->doPrintMemoryInfo = kTRUE;
    StMcParameterDB* parameterDB = StMcParameterDB::instance();
    // TPC
    parameterDB->setXCutTpc(.6); // 6 mm
    parameterDB->setYCutTpc(.6); // 6 mm
    parameterDB->setZCutTpc(.6); // 6 mm
    parameterDB->setReqCommonHitsTpc(10); // Require 10 hits in common for tracks
    // FTPC
    parameterDB->setRCutFtpc(.3); // 3 mm
    parameterDB->setPhiCutFtpc(5*(3.1415927/180.0)); // 5 degrees
    parameterDB->setReqCommonHitsFtpc(3); // Require 3 hits in common for tracks
    // SVT
    parameterDB->setReqCommonHitsSvt(0); // Require 0 hits in common for tracks
    parameterDB->setReqCommonHitsSsd(0); // Require 0 hits in common for tracks

    StEmcAssociationMaker *emcAssociation = new StEmcAssociationMaker(); // EMC association maker
    emcAssociation->setPrint(print);

    StMuDstMaker *muDstMaker = new StMuDstMaker("muDstMaker");

    ///////////////////////////////////////////////////////////////
    //
    // put your analysis maker here
    //
    ///////////////////////////////////////////////////////////////
    
    cout<<" ========================00000000000000000000=============================="<<endl;
    TString outputname = file;
    TString embedrun = file;
    int embedRunIndex = embedrun.Index("_",0);
    embedrun.Remove(0,embedRunIndex+1);
    embedRunIndex = embedrun.Index("P11id",0);
    embedrun.Remove(embedRunIndex);
    cout<<"embedrun:   "<<embedrun<<endl;
    int fileBeginIndex = outputname.Index("st_physics_adc",0);
    outputname.Remove(0,fileBeginIndex+15);
    //fileBeginIndex = outputname.Index("st_physics",0);
    //outputname.Remove(0,fileBeginIndex);
    cout<<" output: "<< outputname<<endl;
    outputname.Prepend("_");
    outputname.Prepend(embedrun);
    outputname.Prepend("emb");
    outputname.ReplaceAll("/","x");
    //outputname.Prepend("st_zerobias");
    outputname.ReplaceAll("event","myminimc");
    outputname.Prepend("/");
    outputname.Prepend(outDir);
    cout<<" Output: "<<outputname<<endl;
    //StEmcTriggerMaker *bemcTrigger = new StEmcTriggerMaker("bemctrigger");
    StTriggerSimuMaker *trigSimuMaker = new StTriggerSimuMaker("StarTrigSimu");
    trigSimuMaker->setMC(1);//0/1 == Using Real/Simulation data files
    trigSimuMaker->useBemc();
    //trigSimuMaker->useEemc(0);
    trigSimuMaker->bemc->setConfig(2);//1: online;  2: offline;  3: expert
    //trigSimuMaker->setPrint(0);
    //StCheckMaker *emcEmbCheck = new StCheckMaker("emcEmbCheck",outputname);

    StMyElectronMaker *electronMaker = new StMyElectronMaker("electronMaker",outputname);
    //end of my analysis maker
    
    chain->Init();
    int iev = 0;
    int istat = 0; 

/*
    int controlVal = 2;
    controlEmcSimulatorMaker_st* simControl = emcSim->getControlSimulator()->GetTable();
    simControl->keyDB[0] = controlVal;
    simControl->keyDB[1] = 0;
    simControl->keyDB[2] = controlVal;
    simControl->keyDB[3] = controlVal;
*/
    // do the event loop    
    while ( istat!=2 && istat!=3 && istat!=4 && iev<=nevents ) {
        chain->Clear();
	cout << "Start to process event number "<<iev <<endl;
        istat = chain->Make();
        //emcAssociation->printMaps();
        //if(iev%10==0) 
	//cout << "Finished processing event number "<<iev <<endl;
        iev++;
    }
    chain->Finish();

}
