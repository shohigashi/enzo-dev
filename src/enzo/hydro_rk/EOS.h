void my_exit(int status);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
inline void EOS(float &p, float &rho, float &e, float &h, float &cs, float &dpdrho,
	       float &dpde, int eostype, int mode)
  /* 
     eostype: 
       0: ideal gas
       1: polytropic EOS
       2: another polytropic EOS
       3: isothermal 
       4: pseudocooling for Wengen 2 test
       5: Wengen 2 the original
       6: minimum pressure similar (but not equal to)
          http://adsabs.harvard.edu/abs/2004ApJ...606...32R
          equation 4
       7: Burkert & Bodenheimer Test. Fedderath et al. (2010) Eqn 30-31
          https://ui.adsabs.harvard.edu/abs/2010ApJ...713..269F

     mode:  
       1: given p and rho, calculate others.
       2: given rho and e, calculate others.
  */
{

  float mu0 = 1.22, mh0 = 1.67e-24, kboltz0 = 1.381e-16;
  float poverrho;

  if (eostype == 0) {
    
    if (mode == 1) {
      poverrho = p / rho;
      e = poverrho / (Gamma - 1);      
    } else if (mode == 2) {
      p = (Gamma - 1) * rho * e;
      poverrho = p / rho;
    }

    dpdrho = poverrho;
    dpde = (Gamma - 1) * rho;
    h = e + poverrho;
    cs = sqrt(Gamma*poverrho);

  }

  if (eostype == 1) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double c_s = EOSSoundSpeed;
    double rho_cr = EOSCriticalDensity;
    //    c_s /= velu;
    rho_cr /= denu;

    cs = c_s*sqrt(1.0 + EOSGamma*pow(rho/rho_cr, EOSGamma-1.0));
    p = rho*c_s*c_s*(1.0 + pow(rho/rho_cr, EOSGamma-1.0));

    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 2) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double c_s = EOSSoundSpeed;
    double rho_cr = 1.0e-12;
    //    c_s /= velu;
    rho_cr /= denu;
    
    if (rho <= 1e-4*rho_cr) {
      cs = c_s*pow(rho/rho_cr, 0.10);
      p = rho*cs*cs;
    } else if(rho <= rho_cr) {
      cs = c_s*pow(rho/rho_cr, -0.20);
      p = rho*cs*cs;
    }else if(rho <= 1e4*rho_cr){ 
      cs = c_s*pow(rho/rho_cr, 0.10);
      p = rho*cs*cs;
    }else if(rho > 1e4*rho_cr){ 
      cs = c_s*pow(rho/rho_cr, 0.666);
      p = rho*cs*cs;
    }
	
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;

  }

  if (eostype == 3) { // straight isothermal
    cs = EOSSoundSpeed;
    p = rho*cs*cs;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 4) { // Wengen 2 test wants pseudocooling
    cs = EOSSoundSpeed;
    // cooling only to 100 should reduce the resolution requirements
    // for the initial tests
    //    cs  = sqrt(1.e-3 + 1./(1.+pow(rho, 1.5)));
    cs *= sqrt(EOSCriticalDensity + 1./(1.+ rho*sqrt(rho)));
    p = rho * cs*cs;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 5) { // Wengen 2 test wants pseudocooling
    // this is the discontinuous one originally suggested
    cs = EOSSoundSpeed;
    // divided by 1000 is the suggested wengen EOS
    // doing to only 100 should reduce the resolution requirements
    // for the initial tests			
    cs = (rho > 1) ?  cs* sqrt(max(1./(rho*sqrt(rho)), 1.e-3)) : cs ;
    p = rho*cs*cs ;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 6) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double rho_cr = EOSCriticalDensity;
    //    c_s /= velu;
    rho_cr /= denu;
    double Pmin = 1.e4*rho_cr/1.67e-24*1.381e-16/(velu*velu);
    Pmin *= (rho/rho_cr)*(rho/rho_cr); // effective gamma = 2 at high density
    if (mode == 1) {
      poverrho = p / rho;
      e = poverrho / (Gamma - 1);      
    } else if (mode == 2) {
      p = (Gamma - 1) * rho * e + Pmin;
      poverrho = p / rho;
    }
    dpdrho = poverrho;
    
    dpde = (Gamma - 1) * rho;
    
    h = e + poverrho;
    // this is somewhat inconsistent as we are using the gamma = 5/3 everywhere
    cs = sqrt(Gamma*poverrho);

  }

  /* Only works in mode 1 and only updates internal energy */
  if (eostype == 7) {
    if (mode == 2) {
        my_exit(EXIT_FAILURE);
    }
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double rho_cr = 1e-15;
    double new_gamma = 0.0;
    double c_s = EOSSoundSpeed;
    rho_cr /= denu;
    
    if(rho/rho_cr <= 0.25) {
      new_gamma = Gamma;
      cs = c_s*pow(rho/rho_cr, new_gamma);
    }
    else if(rho/rho_cr <= 5.0) {
      new_gamma = 1.1;
      cs = c_s*pow(rho/rho_cr, new_gamma);
    }
    else {
      new_gamma = 4.0/3.0;
      cs = c_s*pow(rho/rho_cr, new_gamma);
    }
    
    p = rho*c_s*cs;
    e = p / ((new_gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 8) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
//    double c_s = EOSSoundSpeed;
    double c_s = sqrt(200*1.00*kboltz0/(mu0*mh0)); 
    double rho_cr = EOSCriticalDensity;
    double rho_cr2 = EOSCriticalDensity2;
    double rho_cr3 = EOSCriticalDensity3;
    double c_s0,c_s1, c_s2,c_s3;
    double gamma1 = EOSGamma;
    double gamma2 = EOSGamma2;
    double gamma3 = EOSGamma3;
    c_s /= velu;
    rho_cr /= denu;
    rho_cr2 /= denu;
    rho_cr3 /= denu;

    c_s0 = c_s * sqrt(1.0 + gamma1*pow(rho_cr2/rho_cr,gamma1-1.0))/sqrt(1.0+gamma1);
    c_s1 = c_s * sqrt(1.0 + gamma2*pow(rho_cr2/rho_cr2,gamma2-1.0));
    c_s2 = c_s * sqrt(1.0 + gamma2*pow(rho_cr3/rho_cr2,gamma2-1.0))+(c_s0-c_s1);
    c_s3 = c_s * sqrt(1.0 + gamma3*pow(rho_cr3/rho_cr3,gamma3-1.0));
    if(rho <= rho_cr){
      cs = c_s;
    }else if(rho <= rho_cr2) {
      cs = c_s*sqrt(1.0 + gamma1*pow(rho/rho_cr, gamma1-1.0))/sqrt(1.0+gamma1);
    }else if(rho <= rho_cr3){ 
      cs = c_s*sqrt(1.0 + gamma2*pow(rho/rho_cr2, gamma2-1.0))+(c_s0-c_s1);
    }else{ 
      cs = c_s * sqrt(1.0 + gamma3*pow(rho/rho_cr3,gamma3-1.0))+(c_s2-c_s3);
    }
    p = rho*cs*cs;
    e = p / ((Gamma-1.0)*rho);
//    dpdrho = p/rho;
//    dpde = (Gamma - 1) * rho;
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;

  }
  if (eostype == 9) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
//    double c_s = EOSSoundSpeed;
    double rho_cr = EOSCriticalDensity;
    double gamma1 = EOSGamma;
    double c_s = sqrt(200*gamma1*kboltz0/(mu0*mh0));
    double c_s1, c_s2, p1, p2;
    //c_s1 = sqrt(400*gamma1*kboltz0/(mu0*mh0));
    c_s /= velu;
    //c_s1 /= velu;
    rho_cr /= denu;
    if(rho < rho_cr) {
     cs = c_s;
      //p = pow(cs,2)*pow(rho/rho_cr,1.0);
    }else {
      //cs = c_s*sqrt(1.0 + gamma1*pow(rho/rho_cr, gamma1-1.0))/sqrt(1.0+gamma1);
     cs = c_s*sqrt(pow(rho/rho_cr,gamma1-1.0));
      //cs = c_s;
      //p = pow(cs,2)*pow(rho/rho_cr,1.09);
	  }
    p = rho*cs*cs;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = p/rho;
    dpde = (Gamma - 1) * rho;
//    dpdrho = 1;
//    dpde = 1;
//    h = e + p/rho;
//    dpdrho = pow(cs,2);
//    dpde = (gamma1-1)*rho;
    h = e + dpdrho;

  }

}
