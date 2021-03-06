//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup Common
 * @brief   Definition for rotation matrix
 * @author  Takaya Miyamoto
 * @since   Tue Sep  1 22:25:06 JST 2015
 */
//--------------------------------------------------------------------------

#include <AnalysisHAL.h>

void anaHAL::rot_matrix_init(int rot_matrix[24][4][4], int rot_character[5][24],
                             const int xSIZE, const int ySIZE, const int zSIZE) {
  int R[24][4][4] = {
    { // E
      { 1, 0, 0, 0},
      { 0, 1, 0, 0},
      { 0, 0, 1, 0},
      { 0, 0, 0, 1}
    },
    { // 6C4 (1)
      { 1, 0, 0, 0},
      { 0, 0,-1, xSIZE},
      { 0, 1, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 6C4 (2)
      { 0, 0, 1, 0},
      { 0, 1, 0, 0},
      {-1, 0, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 6C4 (3)
      { 0,-1, 0, xSIZE},
      { 1, 0, 0, 0},
      { 0, 0, 1, 0},
      { 0, 0, 0, 1}
    },
    { // 6C4 (4)
      { 1, 0, 0, 0},
      { 0, 0, 1, 0},
      { 0,-1, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 6C4 (5)
      { 0, 0,-1, xSIZE},
      { 0, 1, 0, 0},
      { 1, 0, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 6C4 (6)
      { 0, 1, 0, 0},
      {-1, 0, 0, xSIZE},
      { 0, 0, 1, 0},
      { 0, 0, 0, 1}
    },
    { // 3C2 (1)
      { 1, 0, 0, 0},
      { 0,-1, 0, xSIZE},
      { 0, 0,-1, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 3C2 (2)
      {-1, 0, 0, xSIZE},
      { 0, 1, 0, 0},
      { 0, 0,-1, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 3C2 (3)
      {-1, 0, 0, xSIZE},
      { 0,-1, 0, xSIZE},
      { 0, 0, 1, 0},
      { 0, 0, 0, 1}
    },
    { // 8C3 (1)
      { 0, 0, 1, 0},
      { 1, 0, 0, 0},
      { 0, 1, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 8C3 (2)
      { 0,-1, 0, xSIZE},
      { 0, 0,-1, xSIZE},
      { 1, 0, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 8C3 (3)
      { 0, 0,-1, xSIZE},
      { 1, 0, 0, 0},
      { 0,-1, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 8C3 (4)
      { 0,-1, 0, xSIZE},
      { 0, 0, 1, 0},
      {-1, 0, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 8C3 (5)
      { 0, 1, 0, 0},
      { 0, 0, 1, 0},
      { 1, 0, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 8C3 (6)
      { 0, 0, 1, 0},
      {-1, 0, 0, xSIZE},
      { 0,-1, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 8C3 (7)
      { 0, 1, 0, 0},
      { 0, 0,-1, xSIZE},
      {-1, 0, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 8C3 (8)
      { 0, 0,-1, xSIZE},
      {-1, 0, 0, xSIZE},
      { 0, 1, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 6C2 (1)
      { 0, 1, 0, 0},
      { 1, 0, 0, 0},
      { 0, 0,-1, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 6C2 (2)
      { 0,-1, 0, xSIZE},
      {-1, 0, 0, xSIZE},
      { 0, 0,-1, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 6C2 (3)
      {-1, 0, 0, xSIZE},
      { 0, 0,-1, xSIZE},
      { 0,-1, 0, xSIZE},
      { 0, 0, 0, 1}
    },
    { // 6C2 (4)
      {-1, 0, 0, xSIZE},
      { 0, 0, 1, 0},
      { 0, 1, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 6C2 (5)
      { 0, 0, 1, 0},
      { 0,-1, 0, xSIZE},
      { 1, 0, 0, 0},
      { 0, 0, 0, 1}
    },
    { // 6C2 (6)
      { 0, 0,-1, xSIZE},
      { 0,-1, 0, xSIZE},
      {-1, 0, 0, xSIZE},
      { 0, 0, 0, 1}
    }
  };
  
  double C[5][5] = {
    // Order: E, 6C4, 3C2, 8C3, 6C2
    { 1, 1, 1, 1, 1}, // A1
    { 1,-1, 1, 1,-1}, // A2
    { 2, 0, 2,-1, 0}, // E
    { 3, 1,-1, 0,-1}, // T1
    { 3,-1,-1, 0, 1}  // T2
  };
  
  for (    int k=0; k<24; k++) 
    for (  int i=0; i< 4; i++) 
      for (int j=0; j< 4; j++) rot_matrix[k][i][j] = R[k][i][j];
  
  for (int i=0; i<5; i++) {
    for (int j= 0; j< 1; j++) rot_character[i][j] = C[i][0]; // E
    for (int j= 1; j< 7; j++) rot_character[i][j] = C[i][1]; // 6C4
    for (int j= 7; j<10; j++) rot_character[i][j] = C[i][2]; // 3C2
    for (int j=10; j<18; j++) rot_character[i][j] = C[i][3]; // 8C3
    for (int j=18; j<24; j++) rot_character[i][j] = C[i][4]; // 6C2
  }
}
