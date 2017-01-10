#!/usr/bin/env python
# vim: ts=4 sw=4 expandtab

import getopt
import os, os.path
import sys, subprocess
import tempfile, math

_install_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
_install_dir = os.environ.get('SENTIEON_INSTALL_DIR', _install_dir)
_paths = { 'gnuplot.sh': 'gnuplot.sh' }

for k,f in _paths.items():
    for d in ('bin', 'libexec'):
        p = os.path.join(os.path.join(_install_dir, d), f)
        if os.path.exists(p):
            _paths[k] = p
            break

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# output to pdf
class GPWriter:
    def __init__(self, out_gp, output, kwargs):
        self.gpf = open(out_gp, 'w')
        self.gpf.write('set terminal pdf size 10in, 7.5in enhanced\n')
        self.gpf.write('set output \"%s\"\n' % output)
        self.kwargs = kwargs
        self.tmpdir = kwargs.get('tmpdir', '/tmp')
        self.tmpfiles = []
        self.enableLogo = kwargs.get('logo', 'true') == 'true'
        self.logoOn = False
    # print log only once per page for multiplot
    # need turn off after the first graph for multiplot
    def printLogo(self):
        if not self.enableLogo:
            return
        self.gpf.write('set label 1 "powered by Sentieon^{/*0.75 TM}" at screen 0.98,0.01 right font "sans,10" tc rgb "#00008B"\n') #darkblue
        self.logoOn = True

    def unprintLogo(self):
        if self.logoOn:
            self.gpf.write('unset label 1\n')
        self.logoOn = False
    def setScreenMargin(self):
        self.gpf.write('set lmargin at screen 0.1\nset rmargin at screen 0.9\nset bmargin at screen 0.1\nset tmargin at screen 0.9\n')
    def close(self):
        self.gpf.close()

    def clean(self):
        for a in self.tmpfiles:
            try:
                os.remove(a)
            except OSError:
                pass

def sharedLabelConvention(xlabel, ylabel, k, n, direction):
    if n == 1:
        return xlabel, ylabel
    xlab = ' '
    ylab = ' '
    if direction == 'x':
        if k == 0:
            ylab = ylabel
        elif k == (n-1)/2:
            xlab = xlabel
    else:
        if k == n-1:
            xlab = ylabel
        elif k == (n-1)/2:
            ylab = ylabel
    return xlab, ylab

class BQSR_GPWriter(GPWriter):
    def __init__(self, out_gp, csv, output, kwargs):
        GPWriter.__init__(self, out_gp, output, kwargs)
        observ_thr = kwargs.get('observations_thr', 10000)
        cutoff = math.floor(math.log10(observ_thr))
        self.expand = int(kwargs.get('expand', '0'))
        self.gpf.write('set datafile separator ","\n')
        self.gpf.write('unset colorbox\n')
        self.gpf.write('unset xtics\nset tic scale 0\nunset border\n')
        self.gpf.write('set palette model RGB defined (-10 "red", -%d "white", %d "white", 10 "blue")\n' % (cutoff, cutoff))
        if self.expand == 2:
            self.printLogo()
            self.setScreenMargin()
        if self.expand == 0:
            self.printLogo()
            self.gpf.write('set multiplot layout 3,3 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')
        self.EventTypes = ['Base Substitution', 'Base Insertion', 'Base Deletion']
        disp_dir = 'x'
        if self.expand == 1:
            disp_dir = 'y'

        CovariateNames = ['QualityScore', 'Cycle', 'Context']
        covIndex = 0
        xlabel = 'Reported Quality Score'
        ylabel = 'Empirical Quality Score'
        fields = [8, 7, 5]
        if self.expand == 1:
            self.printLogo()
            self.gpf.write('set multiplot layout 3,1 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')
        self.plot_scattered_3(csv, CovariateNames[covIndex], xlabel, ylabel, fields, cutoff, 'f(x)=x', 0, 50, 0, 50, disp_dir)

        covIndex = 1
        xlabel = 'Cycle Covariate'
        ylabel = 'Quality Score Accuracy'
        fields = [2, 9, 5]
        if self.expand == 1:
            self.printLogo()
            self.gpf.write('unset multiplot\n')
            self.gpf.write('set multiplot layout 3,1 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')        
        self.plot_scattered_3(csv, CovariateNames[covIndex], xlabel, ylabel, fields, cutoff, 'f(x)=0', -100, 100, -10, 10, disp_dir)
        covIndex = 2
        xlabel = 'Context Covariate'
        ylabel = 'Quality Score Accuracy'
        fields = [0, 9, 5]
        if self.expand == 1:
            self.printLogo()
            self.gpf.write('unset multiplot\n')
            self.gpf.write('set multiplot layout 3,1 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')        
        self.plot_scattered_baselabel_3(csv, CovariateNames[covIndex], xlabel, ylabel, fields, cutoff, 'f(x)=0', 1, 64, -10, 10, disp_dir)

        # print Context cov in another sheet for more detail
        if self.expand == 0:
            self.printLogo()
            self.gpf.write('unset multiplot\n')
            self.gpf.write('set multiplot layout 3,3 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')
        covIndex = 0
        xlabel = 'Quality Score Covariate'
        ylabel = 'Normalized Observations'
        fields = [8, 5]
        if self.expand == 1:
            self.printLogo()
            self.gpf.write('unset multiplot\n')
            self.gpf.write('set multiplot layout 3,1 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')        
        self.plot_hist_3(csv, CovariateNames[covIndex], xlabel, ylabel, fields, 0, '', 0, 55, 0, 0., disp_dir)

        covIndex = 1
        xlabel = 'Cycle Covariate'
        ylabel = 'Mean Quality Score'
        fields = [2, 8, 5]
        if self.expand == 1:
            self.printLogo()
            self.gpf.write('unset multiplot\n')
            self.gpf.write('set multiplot layout 3,1 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')        
        self.plot_scattered_3(csv, CovariateNames[covIndex], xlabel, ylabel, fields, cutoff, '', -100, 100, 0, 50, disp_dir)

        xlabel = 'Context Covariate'
        ylabel = 'Mean Quality Score'
        fields = [0, 7, 5]
        covIndex = 2
        if self.expand == 1:
            self.printLogo()
            self.gpf.write('unset multiplot\n')
            self.gpf.write('set multiplot layout 3,1 title "Covariates Before(red)/After(blue) Recalibration" font "sans,12"\n')        
        self.plot_scattered_baselabel_3(csv, CovariateNames[covIndex], xlabel, ylabel, fields, cutoff, '', 1, 64, 20, 50, disp_dir)

    def plot_scattered_3(self, csv, cov, xlabel, ylabel, fields, cutoff, f_str, xmin, xmax, ymin, ymax, direction = 'x'):
        self.gpf.write('set xrange [%f:%f]\nset yrange [%f:%f]\n' % (xmin, xmax, ymin, ymax))
        self.gpf.write('set xtics font "sans,10"\nset ytics font "sans,10"\n')
        self.gpf.write('%s\n' % f_str)
        k = 0
        narray = 3
        if self.expand == 2:
            narray = 1
        for event in self.EventTypes:
            self.gpf.write('set title "%s"\n' % event)
            xlab, ylab = sharedLabelConvention(xlabel, ylabel, k, narray, direction)
            k += 1
            self.gpf.write('set xlabel "%s"; set ylabel "%s"\n' % (xlab, ylab))            
            src1 = '"<(grep \'%s.*%s.*%s\' %s)"' % (cov, event, 'Before', csv)
            src2 = '"<(grep \'%s.*%s.*%s\' %s)"' % (cov, event, 'After', csv)
            self.gpf.write('plot %s u %d:(int(log10($%d)) >= %d ? $%d : NaN):(-int(log10($5))) notitle w p pt 7 ps .5 palette, \\\n' % (src1, fields[0], fields[2], cutoff, fields[1]))
            self.gpf.write('    %s u %d:(int(log10($%d)) >= %d ? $%d : NaN):(int(log10($5))) notitle w p pt 7 ps .5 palette, \\\n' % (src2, fields[0], fields[2], cutoff, fields[1]))
            self.gpf.write('    f(x) notitle\n')
            if self.expand != 2 and k == 1:
                self.unprintLogo()

    def plot_scattered_baselabel_3(self, csv, cov, xlabel, ylabel, fields, cutoff, f_str, xmin, xmax, ymin, ymax, direction = 'x'):
        self.gpf.write('set xrange [%f:%f]\nset yrange [%f:%f]\n' % (xmin, xmax, ymin, ymax))
        self.gpf.write('set xtics font "sans,10"\nset ytics font "sans,10"\n')
        self.gpf.write('%s\n' % f_str)
        k = 0
        narray = 3
        if self.expand == 2:
            narray = 1
        for event in self.EventTypes:
            self.gpf.write('set title "%s"\n' % event)
            if k == 0:
                self.gpf.write('set xtics rotate by -90 font "sans,10"\n')
                self.gpf.write('set xrange [1:16]\n')
            else:
                if direction == 'x':
                    self.gpf.write('set xtics rotate by -90 font "sans,4"\n')
                else:
                    self.gpf.write('set xtics rotate by -90 font "sans,10"\n')
                self.gpf.write('set xrange [1:64]\n')
            xlab, ylab = sharedLabelConvention(xlabel, ylabel, k, narray, direction)
            k += 1
            self.gpf.write('set xlabel "%s"; set ylabel "%s"\n' % (xlab, ylab))            
            src1 = '"<(grep \'%s.*%s.*%s\' %s)"' % (cov, event, 'Before', csv)
            src2 = '"<(grep \'%s.*%s.*%s\' %s)"' % (cov, event, 'After', csv)
            self.gpf.write('plot %s u %s:(int(log10($%d)) >= %d ? $%d : NaN):(-int(log10($5))):xtic(2) notitle w p pt 7 ps .5 palette, \\\n' % (src1, '(column(0))', fields[2], cutoff, fields[1]))
            self.gpf.write('    %s u %s:(int(log10($%d)) >= %d ? $%d : NaN):(int(log10($5))):xtic(2) notitle w p pt 7 ps .5 palette, \\\n' % (src2, '(column(0))', fields[2], cutoff, fields[1]))
            self.gpf.write('    f(x) notitle\n')
            if self.expand != 2 and k == 1:
                self.unprintLogo()

    # mode 0 by max intensity, mode 1 by area
    def normalizeData(self, csv, filter, seperator, fields, outfile, mode = 0):
        cmd = "grep '%s' %s | awk -F'%s' '{print $%d, $%d}' | sort -k 1 -n > %s" % (filter, csv, seperator, fields[0], fields[1], outfile)
        os.system(cmd)
        norm = 0
        with open(outfile) as fp:
            lines = fp.readlines()
            for line in lines:
                a = line.split()
                val = float(a[1])
                if mode == 0:
                    if val > norm:
                        norm = val
                else:
                    norm += val
        if norm == 0:
            return 0

        maxval= 0
        with open(outfile, "w") as fp:
            for line in lines:
                a = line.split()
                val = float(a[1]) / norm
                if val > maxval:
                    maxval = val
                fp.write('%s,%f\n' % (a[0], val))
        return maxval

    def plot_hist_3(self, csv, cov, xlabel, ylabel, fields, cutoff, f_str, xmin, xmax, ymin, ymax, direction = 'x'):
        self.gpf.write('set xrange [%f:%f]\nset yrange [%f:%f]\n' % (xmin, xmax, ymin, ymax))
        self.gpf.write('set xtics font "sans,10"\nset ytics font "sans,10"\n')
        self.gpf.write('unset yrange\n')
        self.gpf.write('set style line 1 lc rgb "red"\n')
        self.gpf.write('set style line 2 lc rgb "blue"\n')
        self.gpf.write('set style fill transparent solid 1 noborder\n')

        tmpfiles = []
        maxval = 0
        for event in self.EventTypes:
            tmp1 = tempfile.mktemp(prefix='SENTIEON_TMP_',suffix='.before',dir=self.tmpdir)
            tmp2 = tempfile.mktemp(prefix='SENTIEON_TMP_',suffix='.after',dir=self.tmpdir)
            val = self.normalizeData(csv, '%s.*%s.*%s'%(cov, event, 'Before'), ',', fields, tmp1)
            if val > maxval:
                maxval = val
            val = self.normalizeData(csv, '%s.*%s.*%s'%(cov, event, 'After'), ',', fields, tmp2)
            if val > maxval:
                maxval = val
            tmpfiles.append(tmp1)
            tmpfiles.append(tmp2)
        if ymax <= 0:
            ymax = maxval
            self.gpf.write('set yrange [%f:%f]\n' % (ymin, ymax))
        k = 0
        narray = 3
        if self.expand == 2:
            narray = 1
        for event in self.EventTypes:
            self.gpf.write('set title "%s"\n' % event)
            xlab, ylab = sharedLabelConvention(xlabel, ylabel, k, narray, direction)
            self.gpf.write('set xlabel "%s"; set ylabel "%s"\n' % (xlab, ylab))                    
            self.gpf.write('plot "%s" u 1:2 w boxes notitle ls 1, \\\n' % (tmpfiles[k*2]))
            self.gpf.write('"%s" u 1:2 w boxes notitle ls 2\n' % (tmpfiles[k*2+1]))
            k += 1
            if self.expand != 2 and k == 1:
                self.unprintLogo()

        self.tmpfiles = self.tmpfiles + tmpfiles

def plot_bqsr(csv, output, kwargs):
    tmpdir = kwargs.get('tmpdir', '/tmp')
    out_gp = tempfile.mktemp(prefix='SENTIEON_TMP_',suffix='.gp',dir=tmpdir)
    writer = BQSR_GPWriter(out_gp, csv, output, kwargs)
    writer.close()
    os.system('sh %s %s' % (_paths['gnuplot.sh'], out_gp))
    if int(kwargs.get('debug','0')) == 0:
        writer.clean()
        os.remove(out_gp)

metrics_keys = ['gc', 'mq', 'qd', 'isize']

class Metrics_GPWriter(GPWriter):
    def __init__(self, out_gp, inputs, output, kwargs):
        GPWriter.__init__(self, out_gp, output, kwargs)
        self.printLogo()
        for k in metrics_keys:
            if k not in inputs:
                continue
            v = inputs[k]
            if k == 'gc':
                self.plot_gc(v)
            elif k == 'mq':
                self.plot_mq(v)
            elif k == 'qd':
                self.plot_qd(v)
            elif k == 'isize':
                self.plot_isize(v)

    def reset(self):
        self.gpf.write("reset\n")
        self.setScreenMargin()
        self.printLogo()

    def plot_gc(self, v):
        window_max = 0
        key = 'GC'
        with open(v) as f:
            found = False
            for line in f:
                if not found:
                    if line.startswith(key):
                        found = True
                else:
                    words = line.split()
                    if len(words) < 2:
                        continue
                    window = int(words[1])
                    if window > window_max:
                        window_max = window
        if window_max <= 0:
            return
        self.reset()
        self.gpf.write("set title 'GC Bias Plot'\n")
        self.gpf.write("set xlabel 'GC% of 100 base windows'\n")
        self.gpf.write("set xrange [0:100]\n")
        self.gpf.write("set ylabel 'Fraction of normalized coverage'\n")
        self.gpf.write("set yrange [0:2]\n")
        self.gpf.write("set y2label 'Mean base quality'\n")
        self.gpf.write("set y2range [0:40]\n")
        self.gpf.write("set y2tics\n")
        self.gpf.write("set key left\n")
        self.gpf.write("set boxwidth 0.5\n")
        self.gpf.write("set style fill solid 1.0\n")
        self.gpf.write("plot '%s' every ::1 \\\n" % v)
        self.gpf.write("using 1:5:6 title 'Normalized Coverage' with errorbars pt 19 ps 1, \\\n")

        self.gpf.write("'' using ($2/%d/2) title 'Windows at GC%%' with boxes, \\\n" % window_max)
        self.gpf.write("'' using 4 axes x1y2 title 'Base Quality at GC%' with lines\n")

    def plot_mq(self, v):
        self.reset()
        self.gpf.write("set title 'Quality by Cycle'\n")
        self.gpf.write("set xlabel 'Cycle'\n")
        self.gpf.write("set ylabel 'Mean Quality'\n")
        self.gpf.write("set yrange [0:]\n")
        self.gpf.write("plot '%s' every ::1 title 'Mean Quality' with impulses\n" % v)

    def plot_qd(self, v):
        self.reset()
        self.gpf.write("set title 'Quality Score Distribution'\n")
        self.gpf.write("set xlabel 'Quality Score'\n")
        self.gpf.write("set ylabel 'Observations'\n")
        self.gpf.write("set boxwidth 0.5\n")
        self.gpf.write("set style fill solid 1.0\n")
        self.gpf.write("plot '%s' every ::1 title 'Quality Scores' with boxes\n" % v)

    def plot_isize(self, v):
        self.reset()
        self.gpf.write("set title 'Insert Size Histogram'\n")
        self.gpf.write("set xlabel 'Insert Size'\n")
        self.gpf.write("set ylabel 'Count'\n")
        self.gpf.write("set key off\n")
        self.gpf.write("plot '%s' every ::1:1 with lines\n" % v)

def plot_metrics(inputs, output, kwargs):
    tmpdir = kwargs.get('tmpdir', '/tmp')
    out_gp = tempfile.mktemp(prefix='SENTIEON_TMP_',suffix='.gp',dir=tmpdir)
    writer = Metrics_GPWriter(out_gp, inputs, output, kwargs)
    writer.close()
    os.system('sh %s %s' % (_paths['gnuplot.sh'], out_gp))
    if int(kwargs.get('debug','0')) == 0:
        writer.clean()
        os.remove(out_gp)

class VQSR_GPWriter(GPWriter):
    def __init__(self, out_gp, inputs, tranchesFile, output, kwargs):
        GPWriter.__init__(self, out_gp, output, kwargs)
        self.gpf.write('set datafile separator ","\n')
        self.ps = float(kwargs.get('point_size', 0.5))
        self.target_titv = float(kwargs.get('target_titv', 2.15))
        self.min_fp_rate = float(kwargs.get('min_fp_rate', 0.001))
        self.printLogo()
        for k in range(0, len(inputs), 2):
            self.reset()
            self.plot_annotation_pair(inputs[k], inputs[k+1])
        if tranchesFile:
            self.reset()
            self.gpf.write('unset multiplot\n')
            self.plot_tranches(tranchesFile)

    def reset(self):
        self.gpf.write("reset\n")
        self.printLogo()
    def plot_subgraph(self, dataf, field, what, neg, pos):
        self.gpf.write('set key left top reverse box Left title "%s"\n' % what)
        self.gpf.write('plot "%s" u 1:2:($%d<0 ? $%d : NaN) title "%s" w p pt 6 ps %f lc rgb "red", "" u 1:2:($%d>0 ? $%d : NaN) title "%s" w p pt 6 ps %f lc rgb "blue"\n'
                       % (os.path.abspath(dataf), field, field, neg, self.ps, field, field, pos, self.ps))

    def plot_annotation_pair(self, dataf, simuf):
        self.gpf.write('set datafile separator ","\n')
        with open(dataf) as f:
            xlabel = None
            ylabel = None
            for line in f:
                words = line.split(',')
                if len(words) <= 2:
                    print("Cannot not get annotation names from file %s" % dataf)
                    return
                xlabel = words[0]
                ylabel = words[1]
                break
        self.gpf.write('set multiplot layout 2, 2 font "sans,12"\n')
        self.gpf.write('set xlabel "%s" noenhanced; set ylabel "%s" noenhanced\n' % (xlabel, ylabel))
        self.gpf.write('unset key\n')
        self.gpf.write('set pm3d\n')
        self.gpf.write('unset border\n')
        self.gpf.write('set cblabel "lod"\n')
        self.gpf.write('set palette model RGB defined (-4 "red",4 "green")\n')
        self.gpf.write('set cbrange [-4:4]\n')
        self.gpf.write('set title "model PDF"\n')
        self.gpf.write('plot "%s"  u 1:2:3 w p pt 5 ps 1 palette lw 10\n'
                       % os.path.abspath(simuf))
        self.unprintLogo()
        self.gpf.write('set title " "\n')
        self.gpf.write('unset colorbox\n')
        self.gpf.write('unset key\n')
        self.plot_subgraph(dataf, 3, 'outcome','filt','kept')
        self.plot_subgraph(dataf, 4, 'training', 'neg', 'pos')
        self.plot_subgraph(dataf, 5,'novelty', 'novel', 'known')

    def plot_tranches(self, tranchesFile):
        self.setScreenMargin()
        self.plot_tranches_titv(tranchesFile)
        self.reset()
        self.plot_tranches_tpfp(tranchesFile)
    # titv vs truth sensitivity
    def plot_tranches_titv(self, tranchesFile):
        self.gpf.write('set datafile separator ","\n')
        self.gpf.write('unset key\n')
        self.gpf.write('set title "Novel Ti/Tv vs truth sensitivy"\n')
        self.gpf.write('set xlabel "%s"; set ylabel "%s"\n' % ('Tranche truth sensitivy', 'Novel Ti/Tv ratio'))
        self.gpf.write('plot "%s" u 11:5 with lp pt 6 lt rgb "blue"\n' % tranchesFile)
    # true positive/false positive vs titv
    def plot_tranches_tpfp(self, tranchesFile):
        tranches = []
        fields = []
        with open(tranchesFile) as f:
            for line in f:
                if line.startswith('#'):
                    continue;
                line = line.rstrip()
                line = line.lstrip()
                words = line.split(',')
                if len(fields) == 0:
                    fields = words
                    continue
                if len(words) != len(fields):
                    continue
                t = {}
                for i in range(0, len(fields)):
                    if is_number(words[i]):
                        t[fields[i]] = float(words[i])
                    else:
                        t[fields[i]] = words[i]
                fp_rate = self.get_titv_fp(t['novelTiTv'], self.target_titv)
                tp_rate = 1 - fp_rate
                t['numTP'] = int(tp_rate * t['numNovel'])
                t['numFP'] = int(t['numNovel'] - t['numTP'])
                tranches.append(t)
        tmpfile = tempfile.mktemp(prefix='SENTIEON_TMP_',suffix='.tranches.dat',dir=self.tmpdir)
        prev_t = None
        with open(tmpfile, 'w') as f:
            f.write('%s/%s,%s,%s,%s,%s\n' % ('Ti/Tv', 'truth', 'cumulativeFP','cumulativeTP','trancheFP','trancheTP'))
            for t in tranches:
                if prev_t:
                    f.write("%.3f/%.1f,%d,%d,%d,%d\n" % (t['novelTiTv'], 100*t['truthSensitivity'], prev_t['numFP'], prev_t['numTP'],
                                                     t['numFP']-prev_t['numFP'], t['numTP']-prev_t['numTP']))
                else:
                    f.write("%.3f/%.1f,%d,%d,%d,%d\n" % (t['novelTiTv'], 100*t['truthSensitivity'],
                                                     0, 0, t['numFP'], t['numTP']))
                prev_t = t

        self.tmpfiles.append(tmpfile)
        self.reset()
        self.setScreenMargin()
        self.gpf.write('set title "%s"\n' % 'tranche FP/TP vs Ti/Tv')
        self.gpf.write('set xlabel "%s"; set ylabel "%s"\n' % ('Ti/Tv ratio/truth sensitivity', 'numNovel'))
        self.gpf.write('set key left top reverse Left box\n')
        self.gpf.write('set datafile separator ","\n')
        self.gpf.write('set style data histogram\n')
        self.gpf.write('set style histogram rowstacked\n')
        self.gpf.write('set style fill solid border rgb "black"\n')
        self.gpf.write('set auto x\n')
        self.gpf.write('plot "%s" using 2:xtic(1) title column(2) lc rgb "red" fillstyle pattern 1, \\\n' % tmpfile)
        self.gpf.write('    "" using 4:xtic(1) title column(4) lc rgb "red", \\\n')
        self.gpf.write('    "" using 5:xtic(1) title column(5) lc rgb "blue", \\\n')
        self.gpf.write('    "" using 3:xtic(1) title column(3) lc rgb "blue" fillstyle pattern 1\n')

    def get_titv_fp(self, titv, target_titv):
        return max(min(1-(titv-0.5)/(target_titv-0.5),1.), self.min_fp_rate)


def plot_vqsr(inputs, tranchesFile, output, kwargs):
    tmpdir = kwargs.get('tmpdir', '/tmp')
    out_gp = tempfile.mktemp(prefix='SENTIEON_TMP_',suffix='.gp',dir=tmpdir)
    tmpfiles = [out_gp]
    # split into multiple files with data and simu files in pairs, assuming separated by '#' lines
    if len(inputs) == 1:
        with open(inputs[0]) as fp:
            count = 0
            tfp = None
            inputs = []
            for line in fp:
                if line.startswith('#'):
                    if tfp:
                        tfp.close()
                        tfp = None
                else:
                    if not tfp:
                        if count == 0:
                            tmpf = tempfile.mktemp(prefix='SENTIEON_TMP_vqsr_data_',suffix='.txt',dir=tmpdir)
                        else:
                            tmpf = tempfile.mktemp(prefix='SENTIEON_TMP_vqsr_simu_',suffix='.txt',dir=tmpdir)
                        tfp = open(tmpf, 'w')
                        count = (count + 1) % 2
                        inputs.append(tmpf)
                    tfp.write(line)
            if tfp:
                 tfp.close()
            tmpfiles = tmpfiles + inputs
    writer = VQSR_GPWriter(out_gp, inputs, tranchesFile, output, kwargs)
    writer.close()
    os.system('sh %s %s' % (_paths['gnuplot.sh'], out_gp))
    if int(kwargs.get('debug','0')) == 0:
        writer.clean()
        for f in tmpfiles:
            os.remove(f)

def usage(argv):
    print ("Syntax: create plot in pdf format")
    print("    %s bqsr -o output.pdf input.csv" % argv[0])
    print("    %s vqsr -o output.pdf plot_file [tranches_file=tranchesFile] [target_titv=2.15] [min_fp_rate=0.001] [point_size=0.5]'" % argv[0])
    print("    %s metrics -o output.pdf gc=gc.txt mq=mq.txt qd=qd.txt isize=insert_size.txt" % argv[0])

def main(argv):
    if len(argv) < 3:
        usage(argv)
        return 0
    fun = argv[1]
    try:
        opts, args = getopt.getopt(argv[2:], 'o:', ['output'])
    except getopt.GetoptError as err:
        print(str(err))
        usage(argv)
        sys.exit(2)
    output = None
    csv = None
    for o, a in opts:
        if o in ('-o', '--output'):
            output = a
        else:
            assert False, 'unhandled option'
    if fun == 'bqsr':
        if (len(args) < 1 or not output):
            usage(argv)
            sys.exit(2)
        csv = args[0]
        kwargs = dict([arg.split('=',1) for arg in args[1:]])
        if not os.path.isfile(csv):
            print("cannot open file %s\n" % csv)
            sys.exit(2)
        plot_bqsr(csv, output, kwargs)
    elif fun == 'metrics':
        kwargs = dict([arg.split('=',1) for arg in args[0:]])
        inputs = {}
        for k in metrics_keys:
            if k in kwargs:
                inputs[k] = kwargs[k]
                if not os.path.isfile(inputs[k]):
                    print("cannot open file %s\n" % inputs[k])
                    sys.exit(2)
        plot_metrics(inputs, output, kwargs)
    elif fun == 'vqsr':
        if (len(args) < 1 or not output):
            usage(argv)
            sys.exit(2)
        k = 0
        inputs = []
        for arg in args:
            if any(c == '=' for c in arg):
                break
            inputs.append(arg)
            k = k + 1
        kwargs = dict([arg.split('=',1) for arg in args[k:]])
        tranchesFile = kwargs.get('tranches_file', None)
        if tranchesFile and not os.path.isfile(tranchesFile):
            print("cannot open file %s\n" % tranchesFile)
            sys.exit(2)

        for inputfile in inputs:
            if not os.path.isfile(inputfile):
                print("cannot open file %s\n" % inputfile)
                sys.exit(2)
        plot_vqsr(inputs, tranchesFile, output, kwargs)
    else:
        print("Unknown funtion")
        usage(argv)
        sys.exit(2)

if __name__ == "__main__":
    sys.stdin.close()
    sys.stdin = open(os.devnull)
    sys.exit(main(sys.argv))
