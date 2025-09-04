FROM ubuntu:latest
FROM python:3.12-slim

# A: install from local
# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Install any needed packages specified in requirements.txt
RUN pip install -e . --no-cache-dir

#####################################
# B: Install Git to clone repository
# RUN apt-get update
# RUN apt-get install -y git

# Currently need to run it in bash as token for gitlab is needed
# Clone the repository from GitLab
# RUN git clone https://github.com/ferdinand-popp/telomerehunter2.git

# Download package, install package, no cache
# RUN pip install telomerehunter_2/. --no-cache-dir

# Clean up package cache
#RUN apt-get clean
#RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
###################################################

# Optional: Run unit tests after installation
# RUN python -m unittest discover -s tests

# Usage make port available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME=TelomereHunter2
ENV PYTHONWARNINGS="ignore::FutureWarning"

# Set the CMD to start a shell session
CMD ["/bin/bash"]