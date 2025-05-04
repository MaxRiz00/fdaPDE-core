for (int c = 0; c < grid.length; ++c) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      triple p1 = grid[c][i][j];
      triple p2 = grid[c][i+1][j];
      triple p3 = grid[c][i+1][j+1];
      triple p4 = grid[c][i][j+1];

      real s1 = scalarGrid[c][i][j];
      real s2 = scalarGrid[c][i+1][j];
      real s3 = scalarGrid[c][i+1][j+1];
      real s4 = scalarGrid[c][i][j+1];
      real meanS = (s1 + s2 + s3 + s4)/4;

      real t = (meanS - minScalar) / (maxScalar - minScalar); // normalize to [0,1]

      //pen colorPen = rgb(t, 0, 1 - t); // gradient: blue (low) to red (high)

      triple[][] quad = new triple[2][2];
      quad[0][0] = p1;
      quad[0][1] = p2;
      quad[1][0] = p4;
      quad[1][1] = p3;

      //surface patch = surface(quad);
      //draw(patch, surfacepen=material(colorPen + opacity(0.7), emissivepen=gray(0.05), specularpen=mediumgray));

      surface patch = surface(quad);

      // Step 1: add to whole merged surface for lighting/shadow
      wholeSurface = surface(wholeSurface, patch);
      pen colorPen = colormap(t);

      // Step 2: also draw the patch with its individual color
      //draw(patch, surfacepen=material(diffusepen =colorPen, emissivepen=colorPen,specularpen=colorPen));
      draw(patch, surfacepen = material(   diffusepen=gray(0.1),
        specularpen=black,
        emissivepen=colorPen));
    }
  }
}