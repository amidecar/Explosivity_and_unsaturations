<p align="center">
![Project Logo](assets/boom.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
explosivity_and_unsaturations
</h1>

<br>


determining if a molecule is explosive, and find it's unsaturations

##  Usage

```python
#if you want to import only one function, in this example balox
from explosiosivity_and_unsaturations import balox

result_1 = balox(data)

#or import the whole package to get acces to every function
import explosivity_and_unsaturations as eau

result_2 = eau.function(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n explosivity_and_unsaturations python=3.10 
```

```
conda activate explosivity_and_unsaturations
(conda_env) $ pip install .
```

If you need jupyter lab, install it 

```
(explosivity_and_unsaturations) $ pip install jupyterlab
```


## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:Chemikarl/explosivity_and_unsaturations`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:Chemikarl/explosivity_and_unsaturations.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(explosivity_and_unsaturations) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



