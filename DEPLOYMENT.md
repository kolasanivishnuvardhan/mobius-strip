# Deploying the Mobius Strip Explorer

This guide explains how to deploy the Mobius Strip Explorer as a web application using Streamlit.

## Local Deployment

To run the application locally:

1. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

2. Run the Streamlit app:
   ```
   streamlit run streamlit_app.py
   ```

3. The app will open in your default web browser at `http://localhost:8501`

## Deploying to Streamlit Cloud (Free)

1. Create a free account on [Streamlit Cloud](https://streamlit.io/cloud)

2. Connect your GitHub account to Streamlit Cloud

3. Create a new GitHub repository and push your code:
   ```
   git init
   git add .
   git commit -m "Initial commit"
   git remote add origin <your-github-repo-url>
   git push -u origin main
   ```

4. In Streamlit Cloud:
   - Click "New app"
   - Select your repository, branch, and the `streamlit_app.py` file
   - Click "Deploy"

5. Your app will be available at a public URL provided by Streamlit Cloud

## Deploying to Other Platforms

### Heroku

1. Create a `Procfile` with:
   ```
   web: streamlit run streamlit_app.py --server.port=$PORT
   ```

2. Deploy to Heroku:
   ```
   heroku create
   git push heroku main
   ```

### Docker

1. Create a `Dockerfile`:
   ```Dockerfile
   FROM python:3.9-slim
   
   WORKDIR /app
   
   COPY . .
   
   RUN pip install -r requirements.txt
   
   EXPOSE 8501
   
   CMD ["streamlit", "run", "streamlit_app.py"]
   ```

2. Build and run the Docker container:
   ```
   docker build -t mobius-strip-app .
   docker run -p 8501:8501 mobius-strip-app
   ```

## Customization

- To change the appearance, modify the Streamlit theme in the `.streamlit/config.toml` file
- To add authentication, consider using Streamlit's built-in authentication or a service like Auth0