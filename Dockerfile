# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# # Expose port 8501 for streamlit
EXPOSE 8501

# Run streamlit when the container launches
CMD ["streamlit", "run", "Help.py"]
