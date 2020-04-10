# Docker configuration file that adds the action script as an entrypoint to
# the Docker image hosted on DockerHub.

# Use the prebuilt image on DockerHub
FROM bwvdnbro/costuum-ci:latest

# Copy the action script into the container
COPY entrypoint.sh /entrypoint.sh

# Make the action script the entrypoint. The container will execute this
# script and then exit.
ENTRYPOINT ["/entrypoint.sh"]
