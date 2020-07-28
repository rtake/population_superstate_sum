# include <bits/stdc++.h>

using namespace std;


typedef struct EQInfo {
  int num;
  double popl; // population

  vector<int> moverz;
  vector< pair<int,double> > mzpair;
  vector< pair<int,double> > ss;

  // vector<int> frg;
  // vector<double> charge; 
} EQInfo;


typedef struct SSInfo {
  double popl;
  map<int,double> eq; // pair<int eqnum, double coeff>
} SSInfo;


double convert_double(string s) {
  int sign, pos = s.find('e');
  double ret = atof(s.substr(0,pos).c_str());

  if(s[pos+1]=='+') { sign = 1; }
  else if(s[pos+1]=='-') { sign = -1; }

  return ret*pow(10,sign*atoi(s.substr(pos+2,2).c_str()));
}



int main(int argc, char *argv[]) {
  int eqcnt=0, ssnum=0;
  double charge_threshold = 0.5;
  char out[256], base[256], simlog[256], eqpopl[256], eqfrg[256], eqcharge[256], line[256], line0[256], line1[256];
  FILE *fp, *fp0, *fp1, *fp_out;

  vector<EQInfo> eqs;
  SSInfo *sss;

  for(int i=0;i<argc;i++) {
    if( strstr(argv[i],"sim.log") ) { sprintf(simlog,"%s",argv[i]); }
    else if( strstr(argv[i],"_EQ_popl.rrm") ) { sprintf(eqpopl,"%s",argv[i]); }
    else if( strstr(argv[i],".frg") ) { sprintf(eqfrg,"%s",argv[i]); }
    else if( strstr(argv[i],".popl") ) { sprintf(eqcharge,"%s",argv[i]); }
  }

/*
  sprintf(simlog,"cys_ch5_rcmconly_50000k_10+16s_sim.log");
  sprintf(eqpopl,"cys_ch5_rcmconly_50000k_10+16s_EQ_popl.rrm","r");
  sprintf(eqfrg,"list.frg");
  sprintf(eqcharge,"list.popl");
*/

  cout << simlog << endl; 
  cout << eqpopl << endl;
  cout << eqfrg << endl;
  cout << eqcharge << endl;

  sprintf(out,"%s.ss_sum",simlog);
  fp_out = fopen(out,"w");

  fp = fopen(eqpopl,"r");
  for(eqcnt=0;;eqcnt++) {
    if( !fgets(line,256,fp) ) { break; }

    EQInfo eq;
    char population[256], *pt;

    pt = strstr(line,":");
    sscanf(pt+2,"%s",population);

    eq.num = eqcnt;
    eq.popl = convert_double( string(population) );

    eqs.push_back(eq);
  }
  fclose(fp);

  fp0 = fopen(eqfrg,"r");
  fp1 = fopen(eqcharge,"r");

  for(int i=0;i<eqcnt;i++) {
    fgets(line0,256,fp0);
    fgets(line1,256,fp1);

    vector<int> m;
    vector<double> z;

    string s;
    stringstream ss0(line0), ss1(line1);
    
    while( getline(ss0,s,',') ) { m.push_back( atoi( s.c_str() ) ); }
    while( getline(ss1,s,',') ) { z.push_back( atof( s.c_str() ) ); }

    for(int j=0;j<(int)m.size()-1;j++) { eqs[i].mzpair.push_back( pair<int,double>(m[j],z[j]) ); }
  }
  fclose(fp0);
  fclose(fp1);

  fp = fopen(simlog,"r");

  while( fgets(line,256,fp) ) {
    if( strstr(line,"Five representative EQs in SS") ) { sscanf(line,"Five representative EQs in SS   %d\n", &ssnum); } // count number of Super-States

    if( strstr(line,"Coefficients of EQs") ) {
      fgets(line,256,fp);
      break;
    }
  }

  sss = new SSInfo [++ssnum];
  for(int i=0;i<ssnum;i++) { sss[i].popl = 0; }

  for(int i=0;i<eqcnt;i++) {
    int sscnt=0;

    while(true) {
      int index=14;

      for(int j=0;j<10;j++, sscnt++, index += 9) {
        if(index >= (int)string(line).size()) { break; }

        double coeff = atof( string(line).substr(index,8).c_str() );
        if(coeff > 0) {
          eqs[i].ss.push_back( pair<int,double>(sscnt,coeff) );
          sss[sscnt].eq[i] = coeff; // load included EQ and coefficient
          sss[sscnt].popl += eqs[i].popl*coeff; // weighted sum of population
        }
      }

      fgets(line,256,fp);
      if( strstr(line,"EQ -") ) { break; }
      if( string(line) == "\n" ) { break; }
    }
  }

  fclose(fp);

  /*  
  for(int i=0;i<eqcnt;i++) {
    for(int j=0;j<(int)eqs[i].ss.size();j++) {
      printf("%d,", eqs[i].ss[j].first);
    }
    printf("\n");
  } // EQ
  */

  /*
  for(int i=0;i<eqcnt;i++) {
    for(int j=0;j<(int)eqs[i].ss.size();j++) {
      printf("%lf,", eqs[i].ss[j].second);
    }
    printf("\n");
  } // coeff
  */

  for(int i=0;i<eqcnt;i++) {
    for(int j=0;j<(int)eqs[i].mzpair.size();j++) {
      pair<int,double> p = eqs[i].mzpair[j];
      if(p.second > charge_threshold) { eqs[i].moverz.push_back(p.first/((int)(p.second+0.5))); }
    }
  }

  /*
  for(int i=0;i<eqcnt;i++) {
    for(int j=0;j<(int)eqs[i].ss.size();j++) {
      if(eqs[i].ss[j].second > 0) {
        for(int k=0;k<(int)eqs[i].moverz.size();k++) { printf("%d ",eqs[i].moverz[k]); }
      } else { printf("0"); }
      printf(",");
    }
    printf("\n");
  } // coeff
  */

  for(int i=0;i<ssnum;i++) {
    set<int> st_moverz; // store m/z included in the SS
    
    for(map<int,double>::iterator itr=sss[i].eq.begin();itr!=sss[i].eq.end();itr++) {
      for(int j=0;j<(int)eqs[itr->first].mzpair.size();j++) {
        pair<int,double> p = eqs[itr->first].mzpair[j];

        if(p.second > charge_threshold) { st_moverz.insert( p.first/((int)(p.second + 0.5)) ); } // add moverz
      }
    } // for each EQ include in the SS, ...

    
    fprintf(fp_out,"SS-%d : %17.12lf, ", i, sss[i].popl);
    for(set<int>::iterator itr=st_moverz.begin();itr!=st_moverz.end();itr++) { fprintf(fp_out, "%d, ", *itr); }
    fprintf(fp_out,"\n");
    
  } 
 
  fclose(fp_out);

  return 0;
}
