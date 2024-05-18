<p align="center">
![Project Logo](assets/boom.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
explosivity_and_unsaturations
</h1>
<h2 align="center">
Determine if a molecules is explosive and find its unsaturations.

<br>

This program is made of two different parts.
The first one can determine if a molecule is explosive, gives the oxygen balance if it is and highlights in the molecule where the explosible are.
The second part calculates the degree if unsaturation and highlights where these unsaturations are in the molecule.

## How to install the package

First, clone the respository on your own device

```
git clone https://github.com/amidecar/Projeej
```

Select the environment in which you want to install the package or,
Create a new environment, you may also give the environment a different name. 

```
conda create -n explosivity_and_unsaturations python=3.10 
```

```
conda activate explosivity_and_unsaturations
cd path/to/your/cloned/repository
(conda_env)  pip install .
```

If you need jupyter lab, install it 

```
(explosivity_and_unsaturations) $ pip install jupyterlab
```

##  Usage

```python
#if you want to import only one function, in this example balox
from explosiosivity_and_unsaturations import balox

result_1 = balox(data)

#or import the whole package to get acces to every function
import explosivity_and_unsaturations as eau

result_2 = eau.function(data)
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



