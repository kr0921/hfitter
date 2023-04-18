

void mygStyleSetup(){
//gStyle//////////////////////////////////
 gStyle->SetOptStat(00000);
 //gStyle->SetOptStat(00000);
 //gStyle->SetOptTitle(0);
 gStyle->SetHistLineWidth(2);
 gStyle->SetLineWidth(3);
 gStyle -> SetFrameLineWidth(3);
 //gStyle->SetLabelSize(0.035,"X");
 gStyle->SetTitleXOffset(1.2);
 gStyle->SetTitleYOffset(1.35);
 gStyle->SetLabelColor(1,"X");
 gStyle->SetLabelColor(1,"Y");
 //gStyle->SetMarkerStyle(20);
 gStyle -> SetMarkerSize(1);
 gStyle -> SetLabelSize(0.045,"X");
 gStyle -> SetLabelSize(0.045,"Y");
 gStyle -> SetTitleSize(0.045,"XY");
 gStyle -> SetPadLeftMargin(0.14);
 gStyle -> SetPadBottomMargin(0.13);
}



TGraph* count_exp(){
  TFile *_file0 = TFile::Open("res2/graph_expected_ldm.v1.root");
  TGraph* all_fgd12_rec0; _file0.GetObject("all_fgd12_rec0",all_fgd12_rec0);

  double ratio_Mxv =1./3.;
  double alpha_d = 0.5;
  double epsilon = 0.001;
  double Y0 = pow(epsilon,2)*alpha_d*pow(ratio_Mxv,4);
  std::cout<<"Y0= "<<Y0<<std::endl;

  //double n_limit = 3;//(15X)  12.;
  //double n_limit = 2.4;//9;
  
  
  double n_limit = 0.;
  //n_limit = 9; // S/sqrt(S+B)=2 stat=2*10^21POT B=12 => S=9
  //n_limit = 2.4; //statX10 S/sqrt(S+B)=2 B=120 => S=24
  //n_limit = 1.65; //statX20 S/sqrt(S+B)=2 B=240 => S=33   =>33./20=1.65

  n_limit = 5.7; // lim 90%: stat=2*10^21POT B=12 => S=5.7
  //n_limit = 1.5; //statX10 lim90% B=120 => S=15
  //n_limit = 1.05; //statX20 lim90% B=240 => S=21   =>21./20=1.05
  
  Double_t x[200], y[200], y1[200];
  int npoints=0;
  float Mv=0.005;
  while(Mv<0.547){
  //while(Mv<0.127){
     float nrec0 = all_fgd12_rec0->Eval(Mv);  
     //calc with no Eff factor
     //float nrec0 = eta_fgd1_exp0->Eval(Mv)+eta_fgd2_exp0->Eval(Mv);
     //if(Mv<0.127)nrec0 += pi0_fgd1_exp0->Eval(Mv)+pi0_fgd2_exp0->Eval(Mv);
     //float Y_limit1 = Y0*n_limit/nrec0;
     float Y_limit = Y0*sqrt(n_limit/nrec0);
     //dump  std::cout<<Mv<<"  "<<Mv/3.<<"  "<<nrec0<<"  "<<Y_limit/*<<"  "<<Y_limit1*/<<std::endl;
     x[npoints] = Mv/3.;
     y[npoints] = Y_limit;
     y1[npoints] = Y_limit;//Y_limit1;
     npoints++;
     Mv += 0.005;
  }
  TGraph* gr = new TGraph(npoints,x,y);
  gr->SetName("YMx_limits");
  gr->SetTitle("YMx_limits");
  return gr;
}





TGraph* read_file(const char file_name[], int irul=0){
  ifstream in;
  in.open(file_name);
    
  Double_t x[200], y[200];
  int npoints=0;
  while (1) {
    Double_t x1,y1;
    in >> x1 >> y1;
    if (!in.good()) break;
    
    if(y1==0)continue;
    x[npoints] = x1/3.;  // LDM mass
    Double_t epsilon = 0.001*pow(y1,1./4.);
    y[npoints]= pow(epsilon,2)*0.5*pow(1./3.,4.);
    if(irul)std::cout<<"  "<<x1/*x[npoints]*/<<"  "<<y[npoints]<<std::endl;
  
    npoints++;
  }  
  
  TGraph* gr = new TGraph(npoints,x,y);
  gr->SetName("YMx_limits");
  gr->SetTitle("YMx_limits");
  return gr;
}


void ana_sensitivity(){
    
  auto gr = count_exp();
  gr->SetLineColor(1);
  gr->SetLineWidth(3);
    
  auto gr0 =  read_file("res_000/results_sys0.dat",1); 
  gr0->SetLineColor(2);
  gr0->SetLineWidth(3);
    
  auto gr1 =  read_file("res_000/results_sys10.dat"); 
  gr1->SetLineColor(4);
  gr1->SetLineWidth(3);
    
  auto gr2 =  read_file("res_000/results_sys20.dat"); 
  gr2->SetLineColor(5);
  gr2->SetLineWidth(3);
    
  auto gr3 =  read_file("res_000/results_sys30v1.dat"); 
  gr3->SetLineColor(7);
  gr3->SetLineWidth(3);
  
  
  mygStyleSetup();  
  //gStyle -> SetTitleSize(0.045,"XY");

  auto c = new TCanvas();
  c->SetLogx();
  c->SetLogy();
  TH1F *hfr = c->DrawFrame(0.001,1e-12,0.6,1e-7);
  hfr->SetYTitle("Y=#epsilon^{2}#alpha_{D}(m_{#chi}/m_{V})^{4}");
  hfr->SetXTitle("M_{#chi}, GeV");
  hfr->SetTitle("ND280 sensitivity to LDM");
  hfr->GetXaxis()->CenterTitle();
  hfr->GetYaxis()->CenterTitle();
  
  //gr->Draw("same"); 
  gr0->Draw("same"); //0% syst
  //gr1->Draw("same"); 
  //gr2->Draw("same"); 
  gr3->Draw("same"); //30% syst
  
  TLegend* leg = new TLegend(0.19,0.59,0.54,0.84);
  leg->SetBorderSize(0);
  //leg->SetHeader("FGD1");
  leg->AddEntry(gr0,"0\% overall syst.","lp");
  leg->AddEntry(gr3,"30\% overall syst.","lp");
  leg->Draw();
  
  c->SaveAs("ana_sensitivity.png");
}

