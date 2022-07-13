double CalcAlpha3D(kx, ky, kz, px, py, pz){ // 3 momentum components for muon (k) and hadronic system (p)
  double ptx = px + kx;
  double pty = py + ky;
  double ptz = CalcPz(kx, ky, kz, px, py, pz);
  double omegax = -kx;
  double omegay = -ky;
  double omegaz = pz - ptz;

  double dot = omegax*ptx + omegay*pty + omegaz*ptz;
  double norm = TMath::Sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz);
  double norm *= TMath::Sqrt(ptx*ptx + pty*pty + ptz*ptz);
  if (norm == 0){
    return -999; // this means that either q or pn are zero and alpha is undefined.
  }
// old python code checking for ill-defined numbers.  Can re-add if needed.
//  if math.fabs(dot/norm)>0.9999999999:
//    if dot < 0:
//      dot += 1e-10
//    if dot > 0:
//      dot -= 1e-10;
  double alpha3D = TMath::acos(dot/norm);
  return alpha3D;
}

double CalcPn(kx, ky, kz, px, py, pz){// 3 momentum components for muon (k) and hadronic system (p)
  double pt = CalcPt(kx, ky, px, py);
  double pzreco = CalcPz(kx, ky, kz, px, py, pz);
  double pn = TMath::Sqrt(pt*pt + pzreco*pzreco);
  return pn;
}

double CalcPt(kx, ky, px, py){
  double ptx = kx + px;
  double pty = ky + py;
  double pt = TMath::Sqrt(ptx*ptx + pty*pty);
  return pt;
}

double CalcPz(kx, ky, kz, px, py, pz){
  double protonMass = 0.938; // in GeV
  double neutronMass = 0.939; 
  double muonMass = 0.1069; // all of these can be made more precise...
  double Eb = 0.0309;
  double MA = 22*neutronMass + 18*protonMass - 0.3481; // initial nucleus mass
  double MAstar = MA - neutronMass + Eb; // final state neutron mass (average binding energy assumed)
  double Epprime = TMath::Sqrt(protonMass*protonMass + px*px + py*py + pz*pz);
  double Ekprime = TMath::Sqrt(muonMass*muonMass + kx*kx + ky*ky + kz*kz);

  double term1 = MA - Epprime - Ekprime + kz + pz;
  double pt = CalcPt(kx, ky, kz, px, py, pz);
  double term2 = (MAstar*MAstar + pt*pt);
  double pzreco = 0.5*(term1 - term2/term1);
  return pzreco;
}

