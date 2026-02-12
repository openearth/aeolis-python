FROM continuumio/miniconda3

# Set the working directory

# WORKDIR /workspace

# Install Python
RUN conda install python -y

# Copy the project files
# COPY . .

# Install the necessary packages

# RUN pip install .

# Clean up
RUN conda clean -afy

# Set the default command
CMD ["bash"]

# print a message to indicate that the container is ready
RUN echo "The container is ready to use but you need to install aeolis manualy. You can run 'pip instal . -e' to run in editable mode."

