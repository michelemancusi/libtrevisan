
int galois_single_multiply(int x, int y, int w);
int galois_split_w8_multiply(int x, int y);
int galois_create_split_w8_tables();
int galois_inverse(int x, int w);
int galois_shift_inverse(int y, int w);
bool *big_galois_mul(bool *a, bool *b, int l, bool *poly_irr);
bool *big_galois_pow(bool *x, int n, int l, bool *poly_irr);
bool one_bit_ext(bool *subseed_content, bool *source, int n, double error, bool **poly_irr);
int galois_power(int x,int y, int t);
void WDcomputeSi(int i, int m, int t, int t_req, size_t *S);
void wd(int m, int t, int t_req, size_t **S);
void BWDcomputeSi(int *ic, int i, int m, int t, int t_req, size_t *S, size_t *Sc);
void bwd(int m, int t, int t_req, size_t **S);
