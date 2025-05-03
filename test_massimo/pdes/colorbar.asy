
// PARTE MIA

settings.render = 4; 
//settings.prc = true; // if true animation 3d active

import three;
import graph3;
import settings;
import plain;



size(400);
//currentprojection = perspective((0,5,-10),up=(0,1,0));
// sphere ,showtarget=true, autoadjust=false, center=true
//currentprojection =orthographic((2,2,2),up=(0,1,0)); 
currentprojection = perspective((15, 10,10),up=(0,0,1)); // sfera 100 toro 010 , up = (1,1,0)
//currentprojection = perspective((5, 5, 5));
//defaultrender = render(merge = true);
//currentlight = Viewport;

// Soft, frontal lighting + light fill from behind

currentlight = light(
  diffuse = new pen[] {gray(1.0), gray(0.6)},     // brighter light
  //specular = new pen[] {gray(0.3), gray(0.2)},    // subtle highlight
  position = new triple[] {(2, 2, 3), (-2, -1, 2)} // same directions
);



// === SETTINGS ===
int num_points_per_curve = 5;
pen interiorEdgePen = gray + 1bp;
pen boundaryEdgePen = blue + 1.2bp;
pen quadPen = lightblue ;

string folder = "./results/" + settings.user;

// eliminate the last caracter of the string
folder = substr(folder, 0, length(folder) - 1);
folder = folder + "/";


// === HELPERS ===
real clamp(real x, real xmin, real xmax) {
  return max(xmin, min(x, xmax));
}

pen colormap(real t) {
  t = clamp(t, 0, 1);

  // Bordeaux and deep blue
  real blueR = 0.15, blueG = 0.2, blueB = 0.7;
  real redR  = 0.5,  redG  = 0.1, redB  = 0.2;
  real grayR = 0.6,  grayG = 0.6, grayB = 0.6;

  real r, g, b;

  if (t < 0.5) {
    real k = t / 0.5;
    r = (1 - k) * blueR + k * grayR;
    g = (1 - k) * blueG + k * grayG;
    b = (1 - k) * blueB + k * grayB;
  } else {
    real k = (t - 0.5) / 0.5;
    r = (1 - k) * grayR + k * redR;
    g = (1 - k) * grayG + k * redG;
    b = (1 - k) * grayB + k * redB;
  }

  return rgb(r, g, b);
}

pen colormap2(real t) {
  t = clamp(t, 0, 1);
  real r, g, b;

  if (t < 0.2) {
    // Blue to Cyan
    real k = t / 0.2;
    r = 0.0;
    g = k;
    b = 1.0;
  }
  else if (t < 0.4) {
    // Cyan to Green
    real k = (t - 0.2) / 0.2;
    r = 0.0;
    g = 1.0;
    b = 1.0 - k;
  }
  else if (t < 0.6) {
    // Green to Yellow
    real k = (t - 0.4) / 0.2;
    r = k;
    g = 1.0;
    b = 0.0;
  }
  else if (t < 0.8) {
    // Yellow to Orange
    real k = (t - 0.6) / 0.2;
    r = 1.0;
    g = 1.0 - 0.5 * k;
    b = 0.0;
  }
  else {
    // Orange to Red
    real k = (t - 0.8) / 0.2;
    r = 1.0;
    g = 0.5 - 0.5 * k;
    b = 0.0;
  }

  return rgb(r, g, b);
}

triple[] loadTriples(string filename) {
  triple[] result;
  file f = input(filename);
  while (!eof(f)) {
    string line = f;
    string[] p = split(line);
    if (p.length >= 3)
      result.push(((real) p[0], (real) p[1], (real) p[2]));
  }
  return result;
}

int[][] loadEdgeList(string filename) {
  int[][] edges;
  file f = input(filename);
  while (!eof(f)) {
    string line = f;
    string[] p = split(line);
    if (p.length >= 2)
      edges.push(new int[] {((int) p[0] ) , ((int) p[1]) }); // Convert to 1-based indexing
  }
  return edges;
}

int[] loadFlags(string filename) {
  int[] flags;
  file f = input(filename);
  while (!eof(f)) {
    string line = f;
    string[] p = split(line);
    if (p.length >= 1)
      flags.push((int) p[0]);
  }
  return flags;
}

// === LOAD DATA ===
triple[] nodes = loadTriples(folder + "nodes.txt");
triple[] nurbs_edges = loadTriples(folder + "edge_refinement.txt");
int[][] edges = loadEdgeList(folder + "edges.txt");
int[] bflags = loadFlags(folder + "boundary_edges.txt");


/*
triple origin = O; // bottom-left corner of the merged surface

real axisLength = 5.0; // adjust as needed

draw(origin -- (origin + (axisLength+2,0,0)), Arrow3(6bp)); label("$x$", origin + (axisLength+2+0.2,0,0),fontsize(20pt));
draw(origin -- (origin + (0,2*axisLength+2,0)), Arrow3(6bp)); label("$y$", origin + (0,2*axisLength+2+0.2,0),fontsize(20pt));
draw(origin -- (origin + (0,0,axisLength/2)), Arrow3(6bp)); label("$z$", origin + (0,0,axisLength/2+0.2),fontsize(20pt));

*/


// === LOAD & PLOT SURFACE PATCHES ===
file surfFile = input(folder + "refined_surface_points.csv");
triple[][][] grid; // [cell][i][j]
real[][][] scalarGrid; // scalar values at each point
int last_cid = -1;
int cid_index = -1;
int N = num_points_per_curve - 1;

real minScalar = 1e9;
real maxScalar = -1e9;

while (!eof(surfFile)) {
  string line = surfFile;
  string[] p = split(line, ",");

  if (p.length >= 7) { // now includes scalar value
    int cid = (int) p[0];
    int i = (int) p[1];
    int j = (int) p[2];
    real x = (real) p[3];
    real y = (real) p[4];
    real z = (real) p[5];
    real s = (real) p[6];

    if (cid != last_cid) {
      grid.push(new triple[N+1][N+1]);
      scalarGrid.push(new real[N+1][N+1]);
      cid_index += 1;
      last_cid = cid;
    }

    grid[cid_index][i][j] = (x, y, z);
    scalarGrid[cid_index][i][j] = s;

    if (s < minScalar) minScalar = 0; //s
    if (s > maxScalar) maxScalar = 1; //s
  }
}

// print the min and max scalar values
write("Min scalar: " + string(minScalar) );
write("Max scalar: " + string(maxScalar) );

  





// === CREATE 2D OVERLAY PICTURE ===
picture colorbar;

real w = 0.2;
real h = 2.0;
pair origin = (2, 1);
int numSteps = 100;

for (int i = 0; i < numSteps; ++i) {
  real t = i / (real)(numSteps - 1);
  pen color = colormap2(t);
  real y0 = t * h;
  real y1 = (t + 1.0/numSteps) * h;

  fill(colorbar, (origin.x, origin.y + y0) -- 
                 (origin.x + w, origin.y + y0) -- 
                 (origin.x + w, origin.y + y1) -- 
                 (origin.x, origin.y + y1) -- cycle, 
       color);
}

real labelOffset = 0.03; // consistent horizontal offset

label(colorbar, scale(2)*format("%g", minScalar), 
      (origin.x + w + labelOffset, origin.y), E);

label(colorbar, scale(2)*format("%g", (minScalar + maxScalar)/2), 
      (origin.x + w + labelOffset, origin.y + h/2), E);

label(colorbar, scale(2)*format("%g", maxScalar), 
      (origin.x + w + labelOffset, origin.y + h), E);

// === ADD 2D PICTURE TO CURRENT OUTPUT ===
add(currentpicture, colorbar, above=true);  // <-- this is the key line




