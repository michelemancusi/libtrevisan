
#define NONE (10)
#define TABLE (11)
#define SHIFT (12)
#define LOGS (13)
#define SPLITW8 (14)



static int prim_poly[33] = 
{ 0, 
/*  1 */     1, 
/*  2 */    07,
/*  3 */    013,
/*  4 */    023,
/*  5 */    045,
/*  6 */    0103,
/*  7 */    0211,
/*  8 */    0435,
/*  9 */    01021,
/* 10 */    02011,
/* 11 */    04005,
/* 12 */    010123,
/* 13 */    020033,
/* 14 */    042103,
/* 15 */    0100003,
/* 16 */    0210013,
/* 17 */    0400011,
/* 18 */    01000201,
/* 19 */    02000047,
/* 20 */    04000011,
/* 21 */    010000005,
/* 22 */    020000003,
/* 23 */    040000041,
/* 24 */    0100000207,
/* 25 */    0200000011,
/* 26 */    0400000107,
/* 27 */    01000000047,
/* 28 */    02000000011,
/* 29 */    04000000005,
/* 30 */    010040000007,
/* 31 */    020000000011, 
/* 32 */    00020000007 };  /* Really 40020000007, but we're omitting the high order bit */

static int mult_type[33] = 
{ NONE, 
/*  1 */   TABLE, 
/*  2 */   TABLE,
/*  3 */   TABLE,
/*  4 */   TABLE,
/*  5 */   TABLE,
/*  6 */   TABLE,
/*  7 */   TABLE,
/*  8 */   TABLE,
/*  9 */   TABLE,
/* 10 */   LOGS,
/* 11 */   LOGS,
/* 12 */   LOGS,
/* 13 */   LOGS,
/* 14 */   LOGS,
/* 15 */   LOGS,
/* 16 */   LOGS,
/* 17 */   LOGS,
/* 18 */   LOGS,
/* 19 */   LOGS,
/* 20 */   LOGS,
/* 21 */   LOGS,
/* 22 */   LOGS,
/* 23 */   SHIFT,
/* 24 */   SHIFT,
/* 25 */   SHIFT,
/* 26 */   SHIFT,
/* 27 */   SHIFT,
/* 28 */   SHIFT,
/* 29 */   SHIFT,
/* 30 */   SHIFT,
/* 31 */   SHIFT,
/* 32 */   SPLITW8 };

static int nw[33] = { 0, (1 << 1), (1 << 2), (1 << 3), (1 << 4), 
(1 << 5), (1 << 6), (1 << 7), (1 << 8), (1 << 9), (1 << 10),
(1 << 11), (1 << 12), (1 << 13), (1 << 14), (1 << 15), (1 << 16),
(1 << 17), (1 << 18), (1 << 19), (1 << 20), (1 << 21), (1 << 22),
(1 << 23), (1 << 24), (1 << 25), (1 << 26), (1 << 27), (1 << 28),
(1 << 29), (1 << 30), (1 << 31), -1 };

static int nwm1[33] = { 0, (1 << 1)-1, (1 << 2)-1, (1 << 3)-1, (1 << 4)-1, 
(1 << 5)-1, (1 << 6)-1, (1 << 7)-1, (1 << 8)-1, (1 << 9)-1, (1 << 10)-1,
(1 << 11)-1, (1 << 12)-1, (1 << 13)-1, (1 << 14)-1, (1 << 15)-1, (1 << 16)-1,
(1 << 17)-1, (1 << 18)-1, (1 << 19)-1, (1 << 20)-1, (1 << 21)-1, (1 << 22)-1,
(1 << 23)-1, (1 << 24)-1, (1 << 25)-1, (1 << 26)-1, (1 << 27)-1, (1 << 28)-1,
(1 << 29)-1, (1 << 30)-1, 0x7fffffff, 0xffffffff };
   
static int *galois_log_tables[33] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

static int *galois_ilog_tables[33] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

static int *galois_mult_tables[33] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

static int *galois_div_tables[33] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

/* Special case for w = 32 */

static int *galois_split_w8[7] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL };

int galois_create_log_tables(int w)
{
  int j, b;

  if (w > 30) return -1;
  if (galois_log_tables[w] != NULL) return 0;
  galois_log_tables[w] = (int *) malloc(sizeof(int)*nw[w]);
  if (galois_log_tables[w] == NULL) return -1; 
  
  galois_ilog_tables[w] = (int *) malloc(sizeof(int)*nw[w]*3);
  if (galois_ilog_tables[w] == NULL) { 
    free(galois_log_tables[w]);
    galois_log_tables[w] = NULL;
    return -1;
  }
  
  for (j = 0; j < nw[w]; j++) {
    galois_log_tables[w][j] = nwm1[w];
    galois_ilog_tables[w][j] = 0;
  } 
  
  b = 1;
  for (j = 0; j < nwm1[w]; j++) {
    if (galois_log_tables[w][b] != nwm1[w]) {
      fprintf(stderr, "Galois_create_log_tables Error: j=%d, b=%d, B->J[b]=%d, J->B[j]=%d (0%o)\n",
              j, b, galois_log_tables[w][b], galois_ilog_tables[w][j], (b << 1) ^ prim_poly[w]);
      exit(1);
    }
    galois_log_tables[w][b] = j;
    galois_ilog_tables[w][j] = b;
    b = b << 1;
    if (b & nw[w]) b = (b ^ prim_poly[w]) & nwm1[w];
  }
  for (j = 0; j < nwm1[w]; j++) {
    galois_ilog_tables[w][j+nwm1[w]] = galois_ilog_tables[w][j];
    galois_ilog_tables[w][j+nwm1[w]*2] = galois_ilog_tables[w][j];
  } 
  galois_ilog_tables[w] += nwm1[w];
  return 0;
}

int galois_logtable_multiply(int x, int y, int w)
{
  int sum_j;

  if (x == 0 || y == 0) return 0;
  
  sum_j = galois_log_tables[w][x] + galois_log_tables[w][y];
  /* if (sum_j >= nwm1[w]) sum_j -= nwm1[w];    Don't need to do this, 
                                   because we replicate the ilog table twice.  */
  return galois_ilog_tables[w][sum_j];
}

int galois_logtable_divide(int x, int y, int w)
{
  int sum_j;
  int z;

  if (y == 0) return -1;
  if (x == 0) return 0; 
  sum_j = galois_log_tables[w][x] - galois_log_tables[w][y];
  /* if (sum_j < 0) sum_j += nwm1[w];   Don't need to do this, because we replicate the ilog table twice.   */
  z = galois_ilog_tables[w][sum_j];
  return z;
}

int galois_create_mult_tables(int w)
{
  int j, x, y, logx;

  if (w >= 14) return -1;

  if (galois_mult_tables[w] != NULL) return 0;
  galois_mult_tables[w] = (int *) malloc(sizeof(int) * nw[w] * nw[w]);
  if (galois_mult_tables[w] == NULL) return -1;
  
  galois_div_tables[w] = (int *) malloc(sizeof(int) * nw[w] * nw[w]);
  if (galois_div_tables[w] == NULL) {
    free(galois_mult_tables[w]);
    galois_mult_tables[w] = NULL;
    return -1;
  }
  if (galois_log_tables[w] == NULL) {
    if (galois_create_log_tables(w) < 0) {
      free(galois_mult_tables[w]);
      free(galois_div_tables[w]);
      galois_mult_tables[w] = NULL;
      galois_div_tables[w] = NULL;
      return -1;
    }
  }

 /* Set mult/div tables for x = 0 */
  j = 0;
  galois_mult_tables[w][j] = 0;   /* y = 0 */
  galois_div_tables[w][j] = -1;
  j++;
  for (y = 1; y < nw[w]; y++) {   /* y > 0 */
    galois_mult_tables[w][j] = 0;
    galois_div_tables[w][j] = 0;
    j++;
  }
  
  for (x = 1; x < nw[w]; x++) {  /* x > 0 */
    galois_mult_tables[w][j] = 0; /* y = 0 */
    galois_div_tables[w][j] = -1;
    j++;
    logx = galois_log_tables[w][x];
    for (y = 1; y < nw[w]; y++) {  /* y > 0 */
      galois_mult_tables[w][j] = galois_ilog_tables[w][logx+galois_log_tables[w][y]]; 
      galois_div_tables[w][j] = galois_ilog_tables[w][logx-galois_log_tables[w][y]]; 
      j++;
    }
  }
  return 0;
}

int galois_ilog(int value, int w)
{
  if (galois_ilog_tables[w] == NULL) {
    if (galois_create_log_tables(w) < 0) {
      fprintf(stderr, "Error: galois_ilog - w is too big.  Sorry\n");
      exit(1);
    }
  }
  return galois_ilog_tables[w][value];
}

int galois_log(int value, int w)
{
  if (galois_log_tables[w] == NULL) {
    if (galois_create_log_tables(w) < 0) {
      fprintf(stderr, "Error: galois_log - w is too big.  Sorry\n");
      exit(1);
    }
  }
  return galois_log_tables[w][value];
}


int galois_shift_multiply(int x, int y, int w)
{
  int prod;
  int i, j, ind;
  int k;
  int scratch[33];

  prod = 0;
  for (i = 0; i < w; i++) {
    scratch[i] = y;
    if (y & (1 << (w-1))) {
      y = y << 1;
      y = (y ^ prim_poly[w]) & nwm1[w];
    } else {
      y = y << 1;
    }
  }
  for (i = 0; i < w; i++) {
    ind = (1 << i);
    if (ind & x) {
      j = 1;
      for (k = 0; k < w; k++) {
        prod = prod ^ (j & scratch[i]);
        j = (j << 1);
      }
    }
  }
  return prod;
}

int galois_single_multiply(int x, int y, int w)
{
  int sum_j;
  int z;

  if (x == 0 || y == 0) return 0;
  
  if (mult_type[w] == TABLE) {
    if (galois_mult_tables[w] == NULL) {
      if (galois_create_mult_tables(w) < 0) {
        fprintf(stderr, "ERROR -- cannot make multiplication tables for w=%d\n", w);
        exit(1);
      }
    }
    return galois_mult_tables[w][(x<<w)|y];
  } else if (mult_type[w] == LOGS) {
    if (galois_log_tables[w] == NULL) {
      if (galois_create_log_tables(w) < 0) {
        fprintf(stderr, "ERROR -- cannot make log tables for w=%d\n", w);
        exit(1);
      }
    }
    sum_j = galois_log_tables[w][x] + galois_log_tables[w][y];
    z = galois_ilog_tables[w][sum_j];
    return z;
  } else if (mult_type[w] == SPLITW8) {
    if (galois_split_w8[0] == NULL) {
      if (galois_create_split_w8_tables() < 0) {
        fprintf(stderr, "ERROR -- cannot make log split_w8_tables for w=%d\n", w);
        exit(1);
      }
    }
    return galois_split_w8_multiply(x, y);
  } else if (mult_type[w] == SHIFT) {
    return galois_shift_multiply(x, y, w);
  }
  fprintf(stderr, "Galois_single_multiply - no implementation for w=%d\n", w);
  exit(1);
}

int galois_multtable_multiply(int x, int y, int w)
{
  return galois_mult_tables[w][(x<<w)|y];
}

int galois_single_divide(int a, int b, int w)
{
  int sum_j;

  if (mult_type[w] == TABLE) {
    if (galois_div_tables[w] == NULL) {
      if (galois_create_mult_tables(w) < 0) {
        fprintf(stderr, "ERROR -- cannot make multiplication tables for w=%d\n", w);
        exit(1);
      }
    }
    return galois_div_tables[w][(a<<w)|b];
  } else if (mult_type[w] == LOGS) {
    if (b == 0) return -1;
    if (a == 0) return 0;
    if (galois_log_tables[w] == NULL) {
      if (galois_create_log_tables(w) < 0) {
        fprintf(stderr, "ERROR -- cannot make log tables for w=%d\n", w);
        exit(1);
      }
    }
    sum_j = galois_log_tables[w][a] - galois_log_tables[w][b];
    return galois_ilog_tables[w][sum_j];
  } else {
    if (b == 0) return -1;
    if (a == 0) return 0;
    sum_j = galois_inverse(b, w);
    return galois_single_multiply(a, sum_j, w);
  }
  fprintf(stderr, "Galois_single_divide - no implementation for w=%d\n", w);
  exit(1);
}

int galois_shift_divide(int a, int b, int w)
{
  int inverse;

  if (b == 0) return -1;
  if (a == 0) return 0;
  inverse = galois_shift_inverse(b, w);
  return galois_shift_multiply(a, inverse, w);
}

int galois_multtable_divide(int x, int y, int w)
{
  return galois_div_tables[w][(x<<w)|y];
}

void galois_w08_region_multiply(char *region,      /* Region to multiply */
                                  int multby,       /* Number to multiply by */
                                  int nbytes,        /* Number of bytes in region */
                                  char *r2,          /* If r2 != NULL, products go here */
                                  int add)
{
  unsigned char *ur1, *ur2, *cp;
  unsigned char prod;
  int i, srow, j;
  unsigned long l, *lp2;
  unsigned char *lp;
  int sol;

  ur1 = (unsigned char *) region;
  ur2 = (r2 == NULL) ? ur1 : (unsigned char *) r2;

/* This is used to test its performance with respect to just calling galois_single_multiply 
  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) ur2[i] = galois_single_multiply(ur1[i], multby, 8);
  } else {
    for (i = 0; i < nbytes; i++) {
      ur2[i] = (ur2[i]^galois_single_multiply(ur1[i], multby, 8));
    }
  }
 */

  if (galois_mult_tables[8] == NULL) {
    if (galois_create_mult_tables(8) < 0) {
      fprintf(stderr, "galois_08_region_multiply -- couldn't make multiplication tables\n");
      exit(1);
    }
  }
  srow = multby * nw[8];
  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) {
      prod = galois_mult_tables[8][srow+ur1[i]];
      ur2[i] = prod;
    }
  } else {
    sol = sizeof(long);
    lp2 = &l;
    lp = (unsigned char *) lp2;
    for (i = 0; i < nbytes; i += sol) {
      cp = ur2+i;
      lp2 = (unsigned long *) cp;
      for (j = 0; j < sol; j++) {
        prod = galois_mult_tables[8][srow+ur1[i+j]];
        lp[j] = prod;
      }
      *lp2 = (*lp2) ^ l;
    }
  }
  return;
}

void galois_w16_region_multiply(char *region,      /* Region to multiply */
                                  int multby,       /* Number to multiply by */
                                  int nbytes,        /* Number of bytes in region */
                                  char *r2,          /* If r2 != NULL, products go here */
                                  int add)
{
  unsigned short *ur1, *ur2, *cp;
  int prod;
  int i, log1, j, log2;
  unsigned long l, *lp2, *lptop;
  unsigned short *lp;
  int sol;

  ur1 = (unsigned short *) region;
  ur2 = (r2 == NULL) ? ur1 : (unsigned short *) r2;
  nbytes /= 2;


/* This is used to test its performance with respect to just calling galois_single_multiply */
/*
  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) ur2[i] = galois_single_multiply(ur1[i], multby, 16);
  } else {
    for (i = 0; i < nbytes; i++) {
      ur2[i] = (ur2[i]^galois_single_multiply(ur1[i], multby, 16));
    }
  }
  return;
  */

  if (multby == 0) {
    if (!add) {
      lp2 = (unsigned long *) ur2;
      ur2 += nbytes;
      lptop = (unsigned long *) ur2;
      while (lp2 < lptop) { *lp2 = 0; lp2++; }
    }
    return;
  }
    
  if (galois_log_tables[16] == NULL) {
    if (galois_create_log_tables(16) < 0) {
      fprintf(stderr, "galois_16_region_multiply -- couldn't make log tables\n");
      exit(1);
    }
  }
  log1 = galois_log_tables[16][multby];

  if (r2 == NULL || !add) {
    for (i = 0; i < nbytes; i++) {
      if (ur1[i] == 0) {
        ur2[i] = 0;
      } else {
        prod = galois_log_tables[16][ur1[i]] + log1;
        ur2[i] = galois_ilog_tables[16][prod];
      }
    }
  } else {
    sol = sizeof(long)/2;
    lp2 = &l;
    lp = (unsigned short *) lp2;
    for (i = 0; i < nbytes; i += sol) {
      cp = ur2+i;
      lp2 = (unsigned long *) cp;
      for (j = 0; j < sol; j++) {
        if (ur1[i+j] == 0) {
          lp[j] = 0;
        } else {
          log2 = galois_log_tables[16][ur1[i+j]];
          prod = log2 + log1;
          lp[j] = galois_ilog_tables[16][prod];
        }
      }
      *lp2 = (*lp2) ^ l;
    }
  }
  return; 
}

/* This will destroy mat, by the way */

void galois_invert_binary_matrix(int *mat, int *inv, int rows)
{
  int cols, i, j, k;
  int tmp;
 
  cols = rows;

  for (i = 0; i < rows; i++) inv[i] = (1 << i);

  /* First -- convert into upper triangular */

  for (i = 0; i < cols; i++) {

    /* Swap rows if we ave a zero i,i element.  If we can't swap, then the 
       matrix was not invertible */

    if ((mat[i] & (1 << i)) == 0) { 
      for (j = i+1; j < rows && (mat[j] & (1 << i)) == 0; j++) ;
      if (j == rows) {
        fprintf(stderr, "galois_invert_matrix: Matrix not invertible!!\n");
        exit(1);
      }
      tmp = mat[i]; mat[i] = mat[j]; mat[j] = tmp;
      tmp = inv[i]; inv[i] = inv[j]; inv[j] = tmp;
    }
 
    /* Now for each j>i, add A_ji*Ai to Aj */
    for (j = i+1; j != rows; j++) {
      if ((mat[j] & (1 << i)) != 0) {
        mat[j] ^= mat[i]; 
        inv[j] ^= inv[i];
      }
    }
  }

  /* Now the matrix is upper triangular.  Start at the top and multiply down */

  for (i = rows-1; i >= 0; i--) {
    for (j = 0; j < i; j++) {
      if (mat[j] & (1 << i)) {
/*        mat[j] ^= mat[i]; */
        inv[j] ^= inv[i];
      }
    }
  } 
}

int galois_inverse(int y, int w)
{

  if (y == 0) return -1;
  if (mult_type[w] == SHIFT || mult_type[w] == SPLITW8) return galois_shift_inverse(y, w);
  return galois_single_divide(1, y, w);
}

int galois_shift_inverse(int y, int w)
{
  int mat[1024], mat2[32];
  int inv[1024], inv2[32];
  int ind, i, j, k, prod;
 
  for (i = 0; i < w; i++) {
    mat2[i] = y;

    if (y & nw[w-1]) {
      y = y << 1;
      y = (y ^ prim_poly[w]) & nwm1[w];
    } else {
      y = y << 1;
    }
  }

  galois_invert_binary_matrix(mat2, inv2, w);

  return inv2[0]; 
}

int *galois_get_mult_table(int w)
{
  if (galois_mult_tables[w] == NULL) {
    if (galois_create_mult_tables(w)) {
      return NULL;
    }
  }
  return galois_mult_tables[w];
}

int *galois_get_div_table(int w) 
{
  if (galois_mult_tables[w] == NULL) {
    if (galois_create_mult_tables(w)) {
      return NULL;
    }
  }
  return galois_div_tables[w];
}

int *galois_get_log_table(int w)
{
  if (galois_log_tables[w] == NULL) {
    if (galois_create_log_tables(w)) {
      return NULL;
    }
  }
  return galois_log_tables[w];
}

int *galois_get_ilog_table(int w)
{
  if (galois_ilog_tables[w] == NULL) {
    if (galois_create_log_tables(w)) {
      return NULL;
    }
  }
  return galois_ilog_tables[w];
}

void galois_w32_region_multiply(char *region,      /* Region to multiply */
                                  int multby,       /* Number to multiply by */
                                  int nbytes,        /* Number of bytes in region */
                                  char *r2,          /* If r2 != NULL, products go here */
                                  int add)
{
  unsigned int *ur1, *ur2, *cp, *ur2top;
  unsigned long *lp2, *lptop;
  int i, j, a, b, accumulator, i8, j8, k;
  int acache[4];

  ur1 = (unsigned int *) region;
  ur2 = (r2 == NULL) ? ur1 : (unsigned int *) r2;
  nbytes /= sizeof(int);
  ur2top = ur2 + nbytes;

  if (galois_split_w8[0]== NULL) {
    if (galois_create_split_w8_tables(8) < 0) {
      fprintf(stderr, "galois_32_region_multiply -- couldn't make split multiplication tables\n");
      exit(1);
    }
  }

  /* If we're overwriting r2, then we can't do better than just calling split_multiply.
     We'll inline it here to save on the procedure call overhead */

  i8 = 0;
  for (i = 0; i < 4; i++) {
    acache[i] = (((multby >> i8) & 255) << 8);
    i8 += 8;
  }
  if (!add) {
    for (k = 0; k < nbytes; k++) {
      accumulator = 0;
      for (i = 0; i < 4; i++) {
        a = acache[i];
        j8 = 0;
        for (j = 0; j < 4; j++) {
          b = ((ur1[k] >> j8) & 255);
          accumulator ^= galois_split_w8[i+j][a|b];
          j8 += 8;
        }
      }
      ur2[k] = accumulator;
    }
  } else {
    for (k = 0; k < nbytes; k++) {
      accumulator = 0;
      for (i = 0; i < 4; i++) {
        a = acache[i];
        j8 = 0;
        for (j = 0; j < 4; j++) {
          b = ((ur1[k] >> j8) & 255);
          accumulator ^= galois_split_w8[i+j][a|b];
          j8 += 8;
        }
      }
      ur2[k] = (ur2[k] ^ accumulator);
    }
  }
  return;

}

void galois_region_xor(           char *r1,         /* Region 1 */
                                  char *r2,         /* Region 2 */
                                  char *r3,         /* Sum region (r3 = r1 ^ r2) -- can be r1 or r2 */
                                  int nbytes)       /* Number of bytes in region */
{
  long *l1;
  long *l2;
  long *l3;
  long *ltop;
  char *ctop;
  
  ctop = r1 + nbytes;
  ltop = (long *) ctop;
  l1 = (long *) r1;
  l2 = (long *) r2;
  l3 = (long *) r3;
 
  while (l1 < ltop) {
    *l3 = ((*l1)  ^ (*l2));
    l1++;
    l2++;
    l3++;
  }
}

int galois_create_split_w8_tables()
{
  int p1, p2, i, j, p1elt, p2elt, index, ishift, jshift, *table;

  if (galois_split_w8[0] != NULL) return 0;

  if (galois_create_mult_tables(8) < 0) return -1;

  for (i = 0; i < 7; i++) {
    galois_split_w8[i] = (int *) malloc(sizeof(int) * (1 << 16));
    if (galois_split_w8[i] == NULL) {
      for (i--; i >= 0; i--) free(galois_split_w8[i]);
      return -1;
    }
  }

  for (i = 0; i < 4; i += 3) {
    ishift = i * 8;
    for (j = ((i == 0) ? 0 : 1) ; j < 4; j++) {
      jshift = j * 8;
      table = galois_split_w8[i+j];
      index = 0;
      for (p1 = 0; p1 < 256; p1++) {
        p1elt = (p1 << ishift);
        for (p2 = 0; p2 < 256; p2++) {
          p2elt = (p2 << jshift);
          table[index] = galois_shift_multiply(p1elt, p2elt, 32);
          index++;
        }
      }
    }
  }
  return 0;
}

int galois_split_w8_multiply(int x, int y)
{
  int i, j, a, b, accumulator, i8, j8;

  accumulator = 0;
  
  i8 = 0;
  for (i = 0; i < 4; i++) {
    a = (((x >> i8) & 255) << 8);
    j8 = 0;
    for (j = 0; j < 4; j++) {
      b = ((y >> j8) & 255);
      accumulator ^= galois_split_w8[i+j][a|b];
      j8 += 8;
    }
    i8 += 8;
  }
  return accumulator;
}



bool *big_galois_mul(bool *a, bool *b, int l, bool *poly_irr){
	
	bool **Q;
	bool *W;
	bool *T;
	int i=0,j=0, k=0, c=0, d=0,e=0,r=0;
	
	
	Q=malloc(l*sizeof(bool*)); 		
	for(i=0;i<l;i++){
		Q[i]=malloc((2*l-1)*sizeof(bool));	
	}
	
	W=malloc((2*l-1)*sizeof(bool));
	T=malloc(l*sizeof(bool));

	for(i=0;i<l;i++){
		T[i]=0;
	}

	for(i=0;i<l;i++){
		for(j=0;j<(2*l-1);j++){
			Q[i][j]=0;
		}
	}	
	
	for(i=0;i<(2*l-1);i++){
		W[i]=0;
	}
	
	for(i=0;i<l;i++){
		for(j=0;j<l;j++){
			Q[i][2*l-2-j-i]=a[l-1-i]&b[l-1-j];
		}
	}	
	
	for(j=0; j<(2*l-1); j++){
		c=0;
		for(i=0; i<l ;i++){
			c=c^Q[i][j];
		}
		W[j]=c;
	}
	
	d=0;
	e=0;
	r=0;
	
	do{
		d=0;
		e=0;	
		for(j=0; j<(2*l-1-r); j++){
			
			W[j+r]=W[j+r]^poly_irr[j];
						
			if(d==0){  							//DA MIGLIORARE
				if(W[j]==0){
					e++;
				}else{
					d++;
				}
			}else{
				}	
		}
		r=e;
	}while((2*l-1-e)>=(l+1));
	
	
	for(i=0;i<l;i++){				//CICLO DA TOGLIERE SE POSSIBILE, QUELLO CHE FA è PRENDERE L'ARRAY W CHE è LUNGO 2*l-1 E "RIMPICCIOLIRLO" FACENDOLO DIVENRTARE LUNGO l PERCHè L'ARRAY CHE RESTITUISCE LA FUNZIONE DEVE ESSERE DELLA STESSA LUNGHEZZA DEI VETTORI IN INGRESSO POICHè STIAMO IN CAMPO FINITO.
		T[l-i-1]=W[2*l-1-i-1];
	}

	
	for(i=0;i<l;i++){
		free(Q[i]);	
	}
	free(Q);
	free(W);
	
	return T;
}


bool *big_galois_pow(bool *x, int n, int l, bool *poly_irr){
    
    int i=0,k=0, j=0;
    bool *y;
   
   	y=malloc(l*sizeof(bool));
    
    for(i=0;i<(l-1);i++){
		y[i]=0;
	}
	y[l-1]=1;
	
	for(i=0;i<l;i++){
		if(x[i]!=0){
    		j++;
		}
	}
	
	if(n==0){
		return y;
	}else{
		if(n==1){
			return x;
		}else{
		
		do{
			if(n%2==0){
				x=big_galois_mul(x,x,l, poly_irr);
				n=n/2;
			}else{
				y=big_galois_mul(x,y,l, poly_irr);
				x=big_galois_mul(x,x,l, poly_irr);
				n=(n-1)/2;
			}
		}while(n>1);
		}
		
		return big_galois_mul(x,y,l, poly_irr);
	}
	
}


int galois_power(int x,int y, int t){
    int a=0;
    int i=0;
	if(y==0){
		if(x==0){
			return 1;
		}else{
		 	return pow(x,y);
		 }
	}else{
		a=x;
        for (i=0; i<(y-1); i++){
			x=galois_single_multiply(x,a,log2(t));
		
		}
        return x;
    }
}

//Definizione funzione che calcola il singolo sottoseed S_i di lunghezza "t_req"
void WDcomputeSi(int i, int m, int t, int t_req, size_t *S){
    int b=0;
    int j;
    int k;
    int c;
    int *alf;
    int mask;
    int a;
    int Sa;

    
    c=ceil(log2(m)/log2(t_req)-1);
    alf=malloc(c*sizeof(int));
    
    mask=(1<<(int)(log2(t)))-1;

    for (j=0; j<=c; j++){
        alf[j]=(i&(mask<<(j*(int)(log2(t)))))>>j*(int)(log2(t));
    }
                                                   
    for(a=0; a<t_req; a++){
        b=0;
        Sa=0;
		
		for(k=0; k<=c; k++){
			b=b+galois_single_multiply(alf[k],galois_power(a,k,t),log2(t));
        }
		 
		Sa=Sa^b;
		Sa=Sa^(a<<(int)(log2(t)));
		S[a]=Sa;
			
	}    

	free(alf);
}

//Definizione funzione weak design che crea una matrice in cui ci sono "m" sottoseed (S_i) di lunghezza ognuno "t_req"
void wd(int m, int t, int t_req, size_t **S){
	int i;
	for(i=0;i<m;i++){
		WDcomputeSi(i, m, t,t_req,S[i]);
	}
	
}


void BWDcomputeSi(int *ic, int i, int m, int t, int t_req, size_t *S, size_t *Sc){
	double r1=2*M_E;
	int l=MAX(1,ceil((log2(m-r1)-log2(t-r1))/(log2(r1)-log2(r1-1))));
	int j=i%l;
	int k1=i/l;
	int h=0; 
	double n0=((double)m/r1-1);
	int m0=ceil(n0);
	
	if(k1!=*ic){
		*ic=k1;
		WDcomputeSi(*ic, m0, t, t_req, S);
		
		for(h=0; h<t_req; h++){
			//printf("S[%d] vale: %d\n",h,S[h]);
			Sc[h]=S[h];
			//printf("Sc[%d] vale: %d\n",h,Sc[h]);
		}
		
	}else{
		
		for(h=0;h<t_req;h++){
			//printf("g prima %d\n",g);
			S[h]=Sc[h]+j*pow(t,2);
			//printf("g dopo %d\n",g);
			//printf("S[%d] vale %d\n",g,S[g]);
		}
	} 
	//printf("\n");
	
}

void bwd(int m, int t, int t_req, size_t **S){
	int i;
	int ic=-1;
	size_t *Sc;
	Sc=malloc(t_req*sizeof(size_t));
	
	for(i=0; i<m; i++){
		BWDcomputeSi(&ic, i, m, t, t_req, S[i] ,Sc);
	}
	
	free(Sc);
}


bool one_bit_ext(bool *subseed_content, bool *source, int n, double error, bool **poly_irr){
	
	int i=0, j=0;
	int l;
	int s;
	bool **c;
	bool *alpha;
	bool *r;
	bool *a;
	bool b=0;

	l=ceil(log2(n)+2*log2(2/error));
	s=ceil((float)n/(float)l);
		
	a=malloc(l*sizeof(bool));
	r=malloc(l*sizeof(bool));
	alpha=malloc(l*sizeof(bool));
	
	
	c=malloc(s*sizeof(bool*)); 								//VEDERE SE è ALLOCATO BENE COME Q;
	for(i=0;i<s;i++){
		c[i]=malloc(l*sizeof(bool));	
	}
	
	bool *total_source=malloc((s*l)*sizeof(bool));
		
	for(i=0; i<n; i++){
		total_source[i]=source[i];
	}

	for(i=n; i<s*l; i++){
		total_source[i]=0;
	}
	
	for(i=0; i<l; i++){
		r[i]=0;
	}
	
	for(i=0; i<s; i++){
		for(j=0; j<l; j++){
			c[i][j]=total_source[i*l+j];
		}
	}

	for(i=0;i<l;i++){
		alpha[i]=subseed_content[i];
	}

	for(j=1; j<=s; j++){
		a=big_galois_mul(c[j-1], big_galois_pow(alpha, j-1, l, poly_irr[l]), l, poly_irr[l]);
		
		for(i=0; i<l; i++){
			r[i]=r[i]^a[i];	
		}
	}
		
	for(i=0; i<l; i++){
		b=b^(subseed_content[i+l]&r[i]);
	}
	
	for(i=0;i<s;i++){
		free(c[i]);	
	}
	free(c);
	
	free(alpha);
	free(r);
	free(a); 
	free(total_source);
	
	return b;
}
