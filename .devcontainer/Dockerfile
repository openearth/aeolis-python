FROM continuumio/miniconda3

# Set the working directory
WORKDIR /workspace

# Install Python 3.11
RUN conda install python=3.11 -y

# Copy the project files
COPY . .

# Install the necessary packages
RUN pip install .

# Clean up
RUN conda clean -afy

# Set the default command
CMD ["bash"]
