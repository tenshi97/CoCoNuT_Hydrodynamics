#ifdef DOUBLE_PRECISION_EOS
#define EOS_PRECISION double
#else
#define EOS_PRECISION float
#endif


  unsigned long head,tail;
  
  int i,j,k;
  
//struct TABLE { 
int hd_eos_nrho;
int hd_eos_ntt;
int hd_eos_nye;
int hd_eos_timestamp;


EOS_PRECISION *hd_eos_lro_ptr, *hd_eos_ltt_ptr, *hd_eos_ye_ptr, *hd_eos_lpr_ptr;
EOS_PRECISION *hd_eos_led_ptr, *hd_eos_cei_ptr, *hd_eos_cni_ptr,*hd_eos_cpi_ptr;
EOS_PRECISION *hd_eos_cui_ptr, *hd_eos_gai_ptr, *hd_eos_sti_ptr, *hd_eos_xxn_ptr;
EOS_PRECISION *hd_eos_xxp_ptr, *hd_eos_xxa_ptr, *hd_eos_xxh_ptr, *hd_eos_xhz_ptr;
EOS_PRECISION *hd_eos_xha_ptr;
//  };

EOS_PRECISION *hd_eos_lro_lc_ptr, *hd_eos_ltt_lc_ptr, *hd_eos_ye_lc_ptr;  
  
/* global values of the table boundaries */  
EOS_PRECISION den_min, den_max, tem_min, tem_max, ye_min, ye_max;

/* local values of the "chunk"-table boundaries */  
EOS_PRECISION den_min_lc, den_max_lc, tem_min_lc, tem_max_lc, ye_min_lc, ye_max_lc;

int lro_offset, ltt_offset, yei_offset, lpr_offset, led_offset;
int cei_offset, cni_offset, cpi_offset, cui_offset, gai_offset;
int sti_offset, xxn_offset, xxp_offset, xxa_offset, xxh_offset;
int xha_offset, xhz_offset;
  
//struct TABLE hd_eos;
    
