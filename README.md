# py2P Use Guide and Replication Info for Paper Results
This library was used to implement the frameworks and obtain the results in the paper:

A Decentralised Transactive Energy Market Considering Physical System Constraints
https://doi.org/10.1109/TSG.2023.3331714

## Corresponding Framework Implementations
PC2.py implements the transactive PC framework proposed by the above paper, corresponding to results for model "New PC".

SC.py implements the transactive SC framework proposed by the above paper, corresponding to results for model "New SC".

py2Pbase includes the SC and PC implementations for the reference paper used for comparison of results/performance, corresponding to results for models SC[8] and PC[8].

## Results Input Variables 
Results were obtained with the following input parameters:

### 15 bus system
|  model  | trade_scale | deltaP | pen_scale | pen_start | clearing | line_reserve_fix |
| --- | --- | --- | --- | --- | --- | --- |
| New PC | 1e-3 | 1 | 1 | 0 | 0 | 1e-3 |

### 33 bus system    
|  model  | trade_scale | deltaP | pen_scale | pen_start | clearing | line_reserve_fix |
| --- | --- | --- | --- | --- | --- | --- |
| New PC | 1e-3 | 0.1 | 0.1 | 0 | 0 | 1e-3 |

### 33 bus constrained system
|  model  | trade_scale | deltaP | pen_scale | pen_start | clearing | line_reserve_fix |
| --- | --- | --- | --- | --- | --- | --- |
| New PC | 1e-3 | 0.1 | 0.1 | 0 | 0 | 40e-3 |
    
### 141 bus system
|  model  | trade_scale | deltaP | pen_scale | pen_start | clearing | line_reserve_fix |
| --- | --- | --- | --- | --- | --- | --- |
| New PC | 1e-3 | 1 | 1 | 0 | 0 | 1e-3 |

# Additional Notes:
- The implemented code includes debugging options which remain in the primary implementation. In PC2, options pen_start, and clearing are kept for debuggin/testing purposes. 
- Refactoring of code may be completed at some point to improve modularity and quality of implementation for further additions/research.
- Questions or other requests regarding the implemented framework/paper can be sent to colton.d.pankhurst@gmail.com.
