FROM python:3.6

RUN pip install numpy==1.13.3

RUN pip install pandas==0.25.3 

CMD ["python"]



