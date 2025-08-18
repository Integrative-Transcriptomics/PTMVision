# syntax=docker/dockerfile:1
FROM python:3.9-slim
WORKDIR /app
COPY requirements.txt requirements.txt
RUN apt-get update 
RUN apt-get install -y libglib2.0-0 libsm6 libxrender1 libxext6
RUN pip3 install -r requirements.txt
COPY . .
EXPOSE 5001
CMD ["gunicorn","--workers", "4","--timeout","600","--graceful-timeout", "300", "--bind", "0.0.0.0:5001","ptmvision:app"]