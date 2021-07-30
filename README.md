# pigments-from-rrs

## Run in-browser

[![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/cisaacstern/pigments-from-rrs/plotter?filepath=plot_example.ipynb)

Use the link above to open a JupyterLab session in your browser. Changes you make to the repo during the session will not be saved. To save changes made during the session, download any changed files and push them to this repo.

> <table><tr>
>   <td><strong>Note:</strong> JupyterLab link above currently points to this PR branch. It was generated using <a href="https://binder.pangeo.io/">binder.pangeo.io</a> as shown in the image at right. The url will need to be changed to reference <code>alisonpchase/pigments-from-rss</code> on branch <code>main</code> before merge.</td>
>    <td><img src="https://raw.githubusercontent.com/cisaacstern/pigments-from-rrs/images/images/Screen%20Shot%202021-07-30%20at%204.15.26%20PM.png" alt="nbgitpuller example"/></td>
> </tr></table>

<br>

## Local setup

Build and activate the environment:
```
conda env create --file environment.yml \
&& conda activate pigments
```

Install the `pigments` conda enviroment as a Jupyter kernel:
```
python -m ipykernel install --user --name pigments
```

Launch JupyterLab:
```
jupyter lab
```

When opening an notebook in a local JupyterLab session, make sure to select the `pigments` kernel as shown below:

<img src="https://raw.githubusercontent.com/cisaacstern/pigments-from-rrs/images/images/recording.gif" alt="nbgitpuller example fields"/>