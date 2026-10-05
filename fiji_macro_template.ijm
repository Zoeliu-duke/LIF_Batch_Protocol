// =====================================================================================
//  Fiji macro for images processed by process_lif.py
//  (process_lif.py writes processed/fiji_macro.ijm from this template, with the file's
//   channels, colours, display ranges and overlay panels filled in below)
// =====================================================================================
//  HOW TO RUN:  Plugins > Macros > Run...  -> choose processed/fiji_macro.ijm
//    Batch  : choose the 'processed' folder; every *_zstack_raw.tif in each group folder
//             (e.g. processed/10x, processed/63x) is processed. These files hold all
//             z-planes, already stitched for tile scans, calibrated in um.
//             Output: processed/fiji/
//    Active : processes the image open in Fiji (e.g. opened with File > Import >
//             Bio-Formats as a Composite hyperstack, Autoscale OFF).
//  Z handling:
//    Sharpest plane (auto)   : same rule as process_lif.py - focus channels summed,
//                              Gaussian 1.5 px, Laplacian, highest variance wins
//                              (scores printed in the Log window)
//    Single plane (I choose) : type the slice number (Fiji counts from 1)
//    Average / Max projection: merged z-stack (Image > Stacks > Z Project)
//  Then: 32-bit -> Median radius 1 (3x3) -> Gaussian sigma 0.8 px -> fixed display
//  ranges (identical for every sample of a group; never press Auto) -> panel figure.
//  NOTE: the median/Gaussian values below must match DENOISE in process_lif.py.
// =====================================================================================

//@@SETTINGS@@

medianRadius = 1;      // Fiji radius 1 = 3x3 (process_lif.py DENOISE median=3)
gaussSigma = 0.8;      // process_lif.py DENOISE gaussian

zModes = newArray("Sharpest plane (auto)", "Single plane (I choose)", "Average projection", "Max projection");
Dialog.create("Batch image processing");
Dialog.addChoice("Mode", newArray("Batch: processed folder", "Active image"));
Dialog.addChoice("Z handling", zModes);
Dialog.addNumber("Slice to keep, 1-based (only for 'I choose')", 1);
Dialog.addChoice("Group (Active mode only)", groups);
Dialog.show();
mode = Dialog.getChoice();
zMode = Dialog.getChoice();
zPick = Dialog.getNumber();
gActive = Dialog.getChoice();

zi = 0;
for (i = 0; i < zModes.length; i++) if (zMode == zModes[i]) zi = i;
tags = newArray("zAuto", "z" + zPick, "zAVG", "zMAX");
tag = tags[zi];
print("\\Clear");
print("Z handling: " + zMode);

if (startsWith(mode, "Batch")) {
    root = getDirectory("Choose the 'processed' folder");
    out = root + "fiji" + File.separator;
    File.makeDirectory(out);
    setBatchMode(true);
    for (g = 0; g < groups.length; g++) {
        dir = root + groups[g] + File.separator;
        if (!File.exists(dir)) continue;
        files = getFileList(dir);
        for (i = 0; i < files.length; i++) {
            if (!endsWith(files[i], "_zstack_raw.tif")) continue;
            open(dir + files[i]);
            base = replace(files[i], "_zstack_raw.tif", "");
            processImage(g, base, out, zi, zPick, tag, names, luts, focus, bars,
                         single_min, single_max, avg_min, avg_max, max_min, max_max,
                         pairActive, pairTitle, pairRecolor, medianRadius, gaussSigma, barPos);
            close("*");
        }
    }
    setBatchMode(false);
    showMessage("Done", "Saved to " + out + "\nFocus scores are in the Log window.");
} else {
    g = 0;
    for (i = 0; i < groups.length; i++) if (groups[i] == gActive) g = i;
    base = getTitle();
    out = getDirectory("Choose where to save the output");
    processImage(g, base, out, zi, zPick, tag, names, luts, focus, bars,
                 single_min, single_max, avg_min, avg_max, max_min, max_max,
                 pairActive, pairTitle, pairRecolor, medianRadius, gaussSigma, barPos);
}

// =====================================================================================
function processImage(g, base, out, zi, zPick, tag, names, luts, focus, bars,
                      single_min, single_max, avg_min, avg_max, max_min, max_max,
                      pairActive, pairTitle, pairRecolor, medianRadius, gaussSigma, barPos) {
    nch = names.length;
    src0 = getImageID();
    getDimensions(w, h, nc, nz, nt);
    if (nc != nch) exit("This image has " + nc + " channels; the settings list " + nch + ".");

    // ---- 1. z step -------------------------------------------------------------------
    if (zi == 2 && nz > 1) {
        run("Z Project...", "projection=[Average Intensity]");
    } else if (zi == 3 && nz > 1) {
        run("Z Project...", "projection=[Max Intensity]");
    } else {
        z = 1;
        if (zi == 0) z = bestZ(src0, nz, nch, focus, base);
        if (zi == 1) z = minOf(maxOf(zPick, 1), nz);
        selectImage(src0);
        run("Duplicate...", "title=work duplicate slices=" + z);
        print(base + ": kept z-plane " + z + " of " + nz);
    }
    rename("work");

    // ---- 2. display ranges for this group and z option ------------------------------------
    mins = newArray(nch); maxs = newArray(nch);
    for (c = 0; c < nch; c++) {
        k = g * nch + c;
        if (zi == 2)      { mins[c] = avg_min[k];    maxs[c] = avg_max[k]; }
        else if (zi == 3) { mins[c] = max_min[k];    maxs[c] = max_max[k]; }
        else              { mins[c] = single_min[k]; maxs[c] = single_max[k]; }
    }
    bar = bars[g];

    // ---- 3. denoise -----------------------------------------------------------------------
    if (bitDepth() != 32) run("32-bit");
    if (medianRadius > 0) run("Median...", "radius=" + medianRadius + " stack");
    if (gaussSigma > 0) run("Gaussian Blur...", "sigma=" + gaussSigma + " stack");
    if (nch > 1) Stack.setDisplayMode("composite");
    setRanges(nch, luts, "", mins, maxs);
    saveAs("Tiff", out + base + "_" + tag + "_fiji_processed.tif");
    rename("work");
    src = getImageID();

    // ---- 4. panels: top = merge + overlays, bottom = single channels ------------------------
    all = ""; for (c = 0; c < nch; c++) all = all + "1";
    np = pairActive.length;
    ncols = maxOf(1 + np, nch);
    titles = newArray(2 * ncols);
    n = 0;
    titles[n] = toRGB(src, nch, all, "", "Merge", n, luts, mins, maxs, bar, barPos); n++;
    for (p = 0; p < np; p++) { titles[n] = toRGB(src, nch, pairActive[p], pairRecolor[p], pairTitle[p], n, luts, mins, maxs, bar, barPos); n++; }
    for (b = 1 + np; b < ncols; b++) { titles[n] = blank(n, w, h); n++; }
    for (c = 0; c < nch; c++) {
        one = ""; for (k = 0; k < nch; k++) { if (k == c) one = one + "1"; else one = one + "0"; }
        titles[n] = toRGB(src, nch, one, "", names[c], n, luts, mins, maxs, bar, barPos); n++;
    }
    for (b = nch; b < ncols; b++) { titles[n] = blank(n, w, h); n++; }

    selectImage("PNL_00");
    saveAs("PNG", out + base + "_" + tag + "_fiji_merge.png");
    rename("PNL_00");
    run("Images to Stack", "name=panels title=PNL_ use");
    for (i = 1; i <= nSlices; i++) { setSlice(i); setMetadata("Label", titles[i - 1]); }
    run("Make Montage...", "columns=" + ncols + " rows=2 scale=0.5 border=4 font=28 label");
    saveAs("PNG", out + base + "_" + tag + "_fiji_panels.png");
}

// ---- sharpest plane: variance of the Laplacian of the blurred sum of focus channels -----
function bestZ(id, nz, nch, focus, base) {
    fc = split(focus, ",");
    best = 1; bestScore = -1; line = "";
    for (z = 1; z <= nz; z++) {
        selectImage(id);
        run("Duplicate...", "title=fz duplicate slices=" + z);
        run("32-bit");
        if (nch > 1) {
            run("Split Channels");
            acc = "C" + fc[0] + "-fz";
            for (j = 1; j < fc.length; j++) {
                imageCalculator("Add 32-bit", acc, "C" + fc[j] + "-fz");
            }
            selectImage(acc);
        }
        run("Gaussian Blur...", "sigma=1.5");
        run("Convolve...", "text1=[0 1 0\n1 -4 1\n0 1 0]");     // Laplacian
        getStatistics(area, mean, min, max, std);
        score = std * std;
        line = line + "  z" + z + "=" + d2s(score, 1);
        if (score > bestScore) { bestScore = score; best = z; }
        if (nch > 1) { for (c = 1; c <= nch; c++) { if (isOpen("C" + c + "-fz")) close("C" + c + "-fz"); } }
        else close("fz");
    }
    print(base + " focus scores:" + line + "  -> keep z" + best);
    return best;
}

function setRanges(nch, luts, recolor, mins, maxs) {
    for (c = 1; c <= nch; c++) {
        if (nch > 1) Stack.setChannel(c);
        run(luts[c - 1]);
        setMinAndMax(mins[c - 1], maxs[c - 1]);
    }
    if (recolor != "") {                       // e.g. "4:Cyan;2:Yellow"
        items = split(recolor, ";");
        for (i = 0; i < items.length; i++) {
            kv = split(items[i], ":");
            c = parseInt(kv[0]);
            if (nch > 1) Stack.setChannel(c);
            run(kv[1]);
            setMinAndMax(mins[c - 1], maxs[c - 1]);
        }
    }
}

function toRGB(src, nch, active, recolor, label, n, luts, mins, maxs, bar, barPos) {
    selectImage(src);
    setRanges(nch, luts, recolor, mins, maxs);
    if (nch > 1) {
        Stack.setActiveChannels(active);
        run("Stack to RGB");                     // uses only the active channels
    } else {
        run("Duplicate...", "title=tmp"); run("RGB Color");
    }
    t = "PNL_" + IJ.pad(n, 2);
    rename(t);
    // scale bar position from process_lif.py (SCALEBAR_POS); "center" = at a selection in the middle
    getDimensions(bw, bh, bc, bz, bt);
    loc = "Lower Right";
    if (barPos == "lower left") loc = "Lower Left";
    if (barPos == "upper right") loc = "Upper Right";
    if (barPos == "upper left") loc = "Upper Left";
    if (barPos == "center") { makeRectangle(bw / 2 - 1, bh / 2, 2, 2); loc = "At Selection"; }
    run("Scale Bar...", "width=" + bar + " height=8 font=24 color=White background=None location=[" + loc + "] bold");
    run("Select None");
    selectImage(src);
    if (nch > 1) {
        all = ""; for (c = 0; c < nch; c++) all = all + "1";
        Stack.setActiveChannels(all);
        setRanges(nch, luts, "", mins, maxs);
    }
    return label;
}

function blank(n, w, h) {
    newImage("PNL_" + IJ.pad(n, 2), "RGB black", w, h, 1);
    return " ";
}
