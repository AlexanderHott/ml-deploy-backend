FROM python:3.12

WORKDIR /code

COPY ./requirements.in /code/requirements.txt

RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

COPY ./backend /code/backend
COPY ./static /code/static
COPY ./templates /code/templates
RUN mkdir -p /code/results

RUN apt-get update
RUN apt-get install curl unzip -y

RUN curl -L https://brandeis.box.com/shared/static/r16dalx2vodecg93da1gb9micrtmcayt --output checkpoints.zip && unzip checkpoints.zip -d /code/

RUN apt-get install libxrender1 -y

RUN echo $(ls -la /code)

CMD ["fastapi", "run", "backend/main.py", "--host", "0.0.0.0", "--port", "80", "--proxy-headers"]

# If running behind a proxy like Nginx or Traefik add --proxy-headers
# CMD ["fastapi", "run", "app/main.py", "--port", "80", "--proxy-headers"]
