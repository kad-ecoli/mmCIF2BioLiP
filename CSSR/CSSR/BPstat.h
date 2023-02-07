/* parameters for RNA base pairing */
#ifndef BPstat_H
#define BPstat_H 1
const float Pm_mu    = 53.755; const float Pm_sd    =42.773; // P[i-1]-P[i]-P[j]-P[j-1]
const float Pp_mu    = 80.233; const float Pp_sd    =20.466; // P[i+1]-P[i]-P[j]-P[j+1]
const float O5m_mu   = 60.047; const float O5m_sd   =38.731; // O5'[i-1]-O5'[i]-O5'[j]-O5'[j-1]
const float O5p_mu   = 96.645; const float O5p_sd   =20.185; // O5'[i+1]-O5'[i]-O5'[j]-O5'[j+1]
const float C5m_mu   = 58.526; const float C5m_sd   =36.222; // C5'[i-1]-C5'[i]-C5'[j]-C5'[j-1]
const float C5p_mu   =100.129; const float C5p_sd   =23.505; // C5'[i+1]-C5'[i]-C5'[j]-C5'[j+1]
const float C4m_mu   = 64.830; const float C4m_sd   =31.640; // C4'[i-1]-C4'[i]-C4'[j]-C4'[j-1]
const float C4p_mu   =115.546; const float C4p_sd   =23.989; // C4'[i+1]-C4'[i]-C4'[j]-C4'[j+1]
const float C3m_mu   = 80.193; const float C3m_sd   =26.156; // C3'[i-1]-C3'[i]-C3'[j]-C3'[j-1]
const float C3p_mu   =133.291; const float C3p_sd   =26.435; // C3'[i+1]-C3'[i]-C3'[j]-C3'[j+1]
const float C2m_mu   = 87.048; const float C2m_sd   =31.557; // C2'[i-1]-C2'[i]-C2'[j]-C2'[j-1]
const float C2p_mu   =147.550; const float C2p_sd   =29.374; // C2'[i+1]-C2'[i]-C2'[j]-C2'[j+1]
const float C1m_mu   = 74.283; const float C1m_sd   =38.789; // C1'[i-1]-C1'[i]-C1'[j]-C1'[j-1]
const float C1p_mu   =134.728; const float C1p_sd   =27.557; // C1'[i+1]-C1'[i]-C1'[j]-C1'[j+1]
const float O4m_mu   = 62.955; const float O4m_sd   =38.254; // O4'[i-1]-O4'[i]-O4'[j]-O4'[j-1]
const float O4p_mu   =113.863; const float O4p_sd   =24.303; // O4'[i+1]-O4'[i]-O4'[j]-O4'[j+1]
const float O3m_mu   = 80.432; const float O3m_sd   =22.356; // O3'[i-1]-O3'[i]-O3'[j]-O3'[j-1]
const float O3p_mu   =137.697; const float O3p_sd   =28.655; // O3'[i+1]-O3'[i]-O3'[j]-O3'[j+1]
const float P44P_mu  = -0.166; const float P44P_sd  =28.901; // P[i]-C4'[i]-C4'[j]-P[j]
const float C4114C_mu=172.313; const float C4114C_sd=43.508; // C4'[i]-C1'[i]-C1'[j]-C4'[j]
const float PP_mu    = 18.564; const float PP_sd    = 0.896; // P[i]-P[j]
const float O5O5_mu  = 16.808; const float O5O5_sd  = 0.774; // O5'[i]-O5'[j]
const float C5C5_mu  = 17.160; const float C5C5_sd  = 0.513; // C5'[i]-C5'[j]
const float C4C4_mu  = 15.018; const float C4C4_sd  = 0.419; // C4'[i]-C4'[j]
const float C3C3_mu  = 13.676; const float C3C3_sd  = 0.498; // C3'[i]-C3'[j]
const float C2C2_mu  = 11.046; const float C2C2_sd  = 0.515; // C2'[i]-C2'[j]
const float C1C1_mu  = 10.659; const float C1C1_sd  = 0.336; // C1'[i]-C1'[j]
const float O4O4_mu  = 13.232; const float O4O4_sd  = 0.399; // O4'[i]-O4'[j]
const float O3O3_mu  = 15.266; const float O3O3_sd  = 0.641; // O3'[i]-O3'[j]
const float NN_mu    =  8.929; const float NN_sd    = 0.250; // N[i]-N[j]
const float aPm_mu   =109.192; const float aPm_sd   =21.448; // <P[i-1]P[i],P[j+1]P[j]>
const float aPp_mu   =109.326; const float aPp_sd   =20.496; // <P[i+1]P[i],P[j-1]P[j]>
const float aO5m_mu  = 97.545; const float aO5m_sd  =18.629; // <O5'[i-1]O5'[i],O5'[j+1]O5'[j]>
const float aO5p_mu  = 97.816; const float aO5p_sd  =18.227; // <O5'[i+1]O5'[i],O5'[j-1]O5'[j]>
const float aC5m_mu  = 94.309; const float aC5m_sd  =19.724; // <C5'[i-1]C5'[i],C5'[j+1]C5'[j]>
const float aC5p_mu  = 94.640; const float aC5p_sd  =19.843; // <C5'[i+1]C5'[i],C5'[j-1]C5'[j]>
const float aC4m_mu  = 81.590; const float aC4m_sd  =17.368; // <C4'[i-1]C4'[i],C4'[j+1]C4'[j]>
const float aC4p_mu  = 81.352; const float aC4p_sd  =18.082; // <C4'[i+1]C4'[i],C4'[j-1]C4'[j]>
const float aC3m_mu  = 71.338; const float aC3m_sd  =17.754; // <C3'[i-1]C3'[i],C3'[j+1]C3'[j]>
const float aC3p_mu  = 70.590; const float aC3p_sd  =18.100; // <C3'[i+1]C3'[i],C3'[j-1]C3'[j]>
const float aC2m_mu  = 59.411; const float aC2m_sd  =20.249; // <C2'[i-1]C2'[i],C2'[j+1]C2'[j]>
const float aC2p_mu  = 58.352; const float aC2p_sd  =20.141; // <C2'[i+1]C2'[i],C2'[j-1]C2'[j]>
const float aC1m_mu  = 64.542; const float aC1m_sd  =19.817; // <C1'[i-1]C1'[i],C1'[j+1]C1'[j]>
const float aC1p_mu  = 63.616; const float aC1p_sd  =20.319; // <C1'[i+1]C1'[i],C1'[j-1]C1'[j]>
const float aO4m_mu  = 78.198; const float aO4m_sd  =17.830; // <O4'[i-1]O4'[i],O4'[j+1]O4'[j]>
const float aO4p_mu  = 77.672; const float aO4p_sd  =18.859; // <O4'[i+1]O4'[i],O4'[j-1]O4'[j]>
const float aO3m_mu  = 70.536; const float aO3m_sd  =17.393; // <O3'[i-1]O3'[i],O3'[j+1]O3'[j]>
const float aO3p_mu  = 69.919; const float aO3p_sd  =17.786; // <O3'[i+1]O3'[i],O3'[j-1]O3'[j]>
const float aPC_mu   = 59.764; const float aPC_sd   =18.397; // <P[i]C4'[i],P[j]C4'[j]>
const float aCC_mu   =165.213; const float aCC_sd   =13.607; // <C4'[i]C1'[i],C4'[j]C1'[j]>
#endif
