## Read the workspaces 

A script called [readWS.py](readWS.py) allows to read the workspaces from the signal and
background of the actual analysis (also commited here)

## How to play the toys

### Cut and count

A template card for this game is [card_cut_n_count.txt](card_cut_n_count.txt), in which
the strings `%SIG%` and `%BKG%` will be replaced in a script.  
The corresponding script is [toyLimit_cut_n_count.py](toyLimit_cut_n_count.py). Just run it:

``` 
python toyLimit_cut_n_count.py 
```

You will get a plot like this:  
![Cut-n-count results](figs/fig_cut_n_count.png)

As we can see the relationship between the limit and &Sqrt(N_bkg) is indeed
linear. However, we also notice that the line does not cross zero. That means that there
is an offset in the scaling of the limits. We are at the level of small number of events
in the background.  From the backgound `mgg` shape `dN/dmgg = 2 events/GeV`. If we
consider two cases of widths, 2 GeV and 3.2 GeV, we get 4 and 6.4 events
correspondingly. So let's take those numbers and look at how the limits scale:

| N_sig | Limit with N_bkg=4 events | Limit with N_Bkg=6.4 events| (L(4)-L(6.4))/L(4) | compare to &Sqrt;6.4/&Sqrt;4 - 1 |
| - | - | - | - | - |
| 2 | 2.71 | 3.23 | 0.193 | 0.265 |
| 4 | 1.36 | 1.62 | 0.193 | 0.265 |
| 6 | 0.90 | 1.08 | 0.199 | 0.265 |

That is, already with the counting experiment the limit change is smaller than 26%.

### Play with shapes 

Template card for this game is [card_shape_n_roll.txt](card_shape_n_roll.txt). Here we
will be replacing `%SIG%`, `%BKG%` and `%SIGMA%` within
[toyLimit_shape_n_roll.py](toyLimit_shape_n_roll.py) script. In the analysis, Nsig = 4
events, N_bkg = 118 events (total number in 100 &lt; mgg &lt; 180 GeV), sigma = 1.0 (with
a bug) or 1.6 GeV (corrected). Signal shape is Gaussian with mean at 125 GeV and variable
&sigma;, while the background is a simple Exponential function.


The toys generated for that game as well as the background fit used are shown here:
![Toy data](figs/fig_gen_bkg.png)

``` 
python toyLimit_shape_n_roll.py 
```

After running the script, you will get a plot like this, which shows how the limits scale with the width of the gaussian:  
![Shape results](figs/fig_scale_with_sigma.png)

Indeed, the scaling is as &Sqrt;&sigma;, nevertheless the ratio at 1.6 to 1.0 GeV sigmas are as follows:

| N_bkg | N_sig | Limit with &sigma;(sig) = 1 GeV | Limit with &sigma;(sig) = 1.6 GeV| (L(1.0)-L(1.6))/L(1.0) | compare to &Sqrt;1.6/&Sqrt;1 - 1 |
|-|-|-|-|-|-|
| 100 | 2 | 3.10 | 3.70 | 0.194 | 0.265 |
| 120 | 2 | 3.33 | 3.98 | 0.197 | 0.265 |
| 140 | 2 | 3.52 | 4.27 | 0.213 | 0.265 |
| 100 | 4 | 1.55 | 1.85 | 0.194 | 0.265 |
| **120** | **4** | **1.66** | **1.99** | **0.197** | **0.265** |
| 140 | 4 | 1.76 | 2.13 | 0.213 | 0.265 |
| 100 | 6 | 1.04 | 1.24 | 0.196 | 0.265 |
| 120 | 6 | 1.11 | 1.33 | 0.205 | 0.265 |
| 140 | 6 | 1.18 | 1.42 | 0.206 | 0.265 |

(in bold is the case closest to our analysis)

Now, let's do a similar thing and instead of using Gasssian for the signal shape we use Double Sided
Crystal Ball functiio with exactly the same parameters as in the analysis workspaces (before and after the bug fix).
Here are the results:  

<nobr>
<img src="figs/fig_gen_CB_bug.png" width="400"/>
<img src="https://github.com/ResonantHbbHgg/bbggLimits/tree/play-the-toy/toyStudy/figs/fig_gen_CB_bug.png" width="200"/>
</nobr>
![Toy data and CB signal (bugged)](figs/fig_gen_CB_bug.png)
![Toy data and CB signal (fixed)](figs/fig_gen_CB_fix.png)

| N_bkg | N_sig | Limit with bug, &sigma;(eff) = 1 GeV | Limit after fix &sigma;(eff) = 1.6 GeV| (L(bug)-L(fix))/L(bug) | compare to &Sqrt;1.6/&Sqrt;1 - 1 |
|-|-|-|-|-|-|
| 100 | 2 | 3.27 | 3.77 | 0.153 | 0.265 |
| 120 | 2 | 3.52 | 4.05 | 0.151 | 0.265 |
| 140 | 2 | 3.73 | 4.33 | 0.159 | 0.265 |
| 100 | 4 | 1.63 | 1.88 | 0.153 | 0.265 |
| **120** | **4** | **1.76** | **2.02** | **0.151** | **0.265** |
| 140 | 4 | 1.87 | 2.16 | 0.159 | 0.265 |
| 100 | 6 | 1.09 | 1.26 | 0.158 | 0.265 |
| 120 | 6 | 1.17 | 1.36 | 0.161 | 0.265 |
| 140 | 6 | 1.24 | 1.44 | 0.164 | 0.265 |

The 15% difference is what we observe in the analysis as well. 

Note: the width of 1.0 and 1.6 GeV are the _effective sigmas_ of the PDF. In the case of
Crystal Ball PDFs with large tails, a simple scaling with the &sigma;_eff breaks.

 
 
 
 
 
 
 
 


