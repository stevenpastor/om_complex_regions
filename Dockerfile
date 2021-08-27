FROM python:3.6

RUN pip install numpy==1.13.3

RUN pip install pandas==0.25.3 

CMD ["python"]

# docker run --rm --mount type=bind,source=“$(pwd)“/data,target=/data --mount type=bind,source=“$(pwd)“/scripts,target=/scripts -it python:3.6 python scripts/initial_check.py data/


