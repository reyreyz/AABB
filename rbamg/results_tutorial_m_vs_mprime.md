# Tutorial code (build_it_yourself): M (dP=0) vs M' (dP>=1), best of the same knob grid, m=180, tol 1e-10

```

stagnation:
  M (dP=0)     W0=direct  :   91.0 WU,  40 cycles, op.cx 2.86   (theta, dP, pre, post, filt) = (0.25, 0, 'F', '-', 0.001)
  M (dP=0)     W0=onepoint:  133.7 WU,  51 cycles, op.cx 2.16   (theta, dP, pre, post, filt) = (0.25, 0, 'F', 'F', '-')
  M' (dP>=1)   W0=direct  :   19.4 WU,   4 cycles, op.cx 3.73   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M' (dP>=1)   W0=onepoint:   19.3 WU,   4 cycles, op.cx 3.72   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M' (dP>=1)   W0=zero    :   22.4 WU,   5 cycles, op.cx 3.44   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.001)

curved:
  M (dP=0)     W0=direct  :   89.9 WU,  41 cycles, op.cx 2.88   (theta, dP, pre, post, filt) = (0.25, 0, 'F', '-', 0.0001)
  M' (dP>=1)   W0=direct  :   21.6 WU,   4 cycles, op.cx 4.14   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M' (dP>=1)   W0=onepoint:   21.1 WU,   4 cycles, op.cx 4.09   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M' (dP>=1)   W0=zero    :   26.7 WU,   5 cycles, op.cx 4.09   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M (dP=0)     W0=onepoint: no configuration converges

jump:
  M (dP=0)     W0=direct  :   55.7 WU,  30 cycles, op.cx 2.28   (theta, dP, pre, post, filt) = (0.25, 0, 'F', '-', 0.0001)
  M' (dP>=1)   W0=direct  :   12.9 WU,   4 cycles, op.cx 2.61   (theta, dP, pre, post, filt) = (0.25, 2, 'FFF', '-', 0.0001)
  M' (dP>=1)   W0=onepoint:   13.0 WU,   4 cycles, op.cx 2.60   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M' (dP>=1)   W0=zero    :   13.0 WU,   4 cycles, op.cx 2.60   (theta, dP, pre, post, filt) = (0.25, 3, 'FFF', '-', 0.0001)
  M (dP=0)     W0=onepoint: no configuration converges
```
