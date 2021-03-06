ó
:XRc           @   s   d  d l  j Z d  d l Z d  d l Z d  d l Z d   Z d d d d d d  Z	 d e
 d d  Z d   Z e d k r e   n  d S(   iÿÿÿÿNc         C   sJ   |  d |  d  } |  d  | d } | d | d  } | | } | | f S(   Ni   iÿÿÿÿi   (    (   t   xt   yt   dif_xt   x_difft   dif_yt   y_diff(    (    s4   /home/ezequiel/.local/lib/pymodules/eTools/pTools.pyt   DiffData   s
    
t
   Polynomialc      
      sÅ  d3 |    f d  } |  r[ g  |  D] } | d ^ q% } g  |  D] } | d ^ qB } n | rq | rq d GHn  | s t j j t j    } n  d } g  t t |   D] } d j |  ^ q« }	 t j	   }
 |
 j
 |  |
 j d  } | j d  | j d	  | j | | d
 d d d d d xH t |	 | |  D]4 \ } } } t j | d | | f d d4 d d q=Wt j d | d d d d d t j   t j d d t d t \ } \      j d   j d   j d    j d! d" d# d$ d% d&  | j d' d  | j
 |    j | | d
 d( d d d d) | } | } x>t rd* GHd+ GHd, GHt d-  } t j d. |  } | d d/ k rºt |  d k rÀt | d  } n d3 } g  } g  } xá t t |   D]Í } t j | t | |   rå| j | |  | j | |  | t |  d k rgt |  d k rg| | | d0 | q²| | d d	 k r²t |  d k r²| | | d0 | g  } g  } q²qåqåWPqX| d d
 k rÎPqXt | d  } t | d  } t |  d1 k rt | d  } n d3 } | | | d !} | | | d !} g  t |  D]1 \ } } | t | | d  k rsd	 n | ^ qH} | | | d0 | qXWt j d2 | d d d d d t j   d3 S(5   sª   
    Linear regression of square fit of some x,y data.
    type, can be 'Polynomial' or 'Chevbyshev'.
    interv, boolean. Weather it creates a total or partial fit.
    c      	      sl   t  |  | d t d | d | \ } } }   j | | d d d d d d	  j | | d d
 d d	 d S(   s+   
        Add partial fits to axes.
        t   forcet   gradt   typet   colort   rt	   linestylet   -t	   linewidths   1.5t   blueN(   t   FitDatat   Truet   plot(   t   x_to_fitt   y_to_fitR	   t   polit   xfitt   yfitt   dyfit(   t   ax1t   ax2(    s4   /home/ezequiel/.local/lib/pymodules/eTools/pTools.pyt   plot_partial_fit   s
    i    i   s2   Must declare either "x" and "y" or two-element sets   -?(\d+(\.\d*)?|\d*\.d+)s   P{0}io   R   R    t   si2   t
   facecolorst   nonet
   edgecolorst   redt   xyt   xytexti   t
   textcoordss   offset pointst   fit_s   .pngt   formatt   pngt   bbox_inchest   tighti   t   sharext   shareys   Energy [eV]s
   Force [nN]s   Stretching distance [\AA]R   s   --R   t   grayR   t   2t   hspacei   R   s*   Range of point numbers? (comma separated):s&   Enter "r"  to fit all remaining pointss   Enter "s"  to stop fittings   --> s   ,|:R   R	   i   s   energy-force_N(   i    i   (    t   Nonet   ost   patht   basenamet   getcwdt   ranget   lenR&   t   pltt   figuret   suptitlet   add_subplott
   set_ylabelt
   set_xlabelt   scattert   zipt   annotatet   savefigt   showt   subplotsR   t   Falset   axhlinet   subplots_adjustt	   raw_inputt   ret   splitt   intt   matcht   strt   appendt	   enumerate(   t   dataR    R   t   sysNameR   R   t   a_datat   real_pt   it   labelst   f0t   ax0t   labelt   axt   ayt   f2t   left_xt   left_yt   selt   selistt   gradoR   R   t   at   bt   e(    (   R   R   s4   /home/ezequiel/.local/lib/pymodules/eTools/pTools.pyt   FitBlock   s    
 .""!
'"	(&Ac   	      C   s)  t  j |  d |  d d t |   d } | sb t t |   d  d } | d k rb d } qb n  | d k r· t  j j j |  | |  } t  j | d	 d	 d   } t  j |  } n? | d
 k rö t  j j j	 |  | |  } t  j j j
 |  } n  | r| | |  | |  f S| | |  f Sd	 S(   s­   
    Return x and y values of best fit for x y data set by square fit of desired
    type.

    type: at this moment 'Polynomial' or 'Chebyshev' are the posible options
    i    iÿÿÿÿt   numi   i   i   i   R   Nt	   Chebyshev(   t   npt   linspaceR5   RH   t
   polynomialt   polyfitt   poly1dt   polyderRc   t   fitt   deriv(	   R    R   R
   R   R	   t   x_newt   coefst   ffitt   ffitd(    (    s4   /home/ezequiel/.local/lib/pymodules/eTools/pTools.pyR   v   s    *c           C   s   d  S(   N(    (    (    (    s4   /home/ezequiel/.local/lib/pymodules/eTools/pTools.pyt   main   s    t   __main__(   t   matplotlib.pyplott   pyplotR6   t   numpyRd   R0   RF   R   R/   Ra   RB   R   Rp   t   __name__(    (    (    s4   /home/ezequiel/.local/lib/pymodules/eTools/pTools.pyt   <module>   s   	d	