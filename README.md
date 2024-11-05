# Smiles ML Model

![image](https://github.com/user-attachments/assets/061b2870-c731-4742-ad3f-c7b8b3cd3abe)

The technoligies I used were
- fastapi: backend framework
- htmx: minimal frontend library to make the website feel faster
- sqlite: simple database

I went with fastapi because the code to run the model was already in python, and it is a fast (who would have guessed) and simple to use.

I went with HTMX over no front end library because the bare html form UX is bad and requires the server to send more data.

## Performance

On my laptop:
- AMD Ryzen 7 7840U w/ Radeon  780M Graphics
- 16GB DDR4 ram

For the home page:

- Success rate: 100.00%
- Total:        9.4908 secs
- Slowest:      0.1110 secs
- Fastest:      0.0410 secs
- Average:      0.0474 secs
- Requests/sec: 1053.6476

For predictions:

- Success rate: 90.30%
- Total:        69.4551 secs
- Slowest:      15.4782 secs
- Fastest:      0.3467 secs
- Average:      3.6887 secs
- Requests/sec: 14.3978

The success rate was bad, but that can be fixed with a queue like [celery](https://docs.celeryq.dev/en/stable/getting-started/introduction.html).
A queuing system would also allow the server to handle variable load on a less powerful (cheaper) CPU at the cost of instant feedback for the user.

- [200] 679 responses
- [500] 224 responses

## Pricing

The project could run off of a $10/mo server for small usages, but it really depends on how many requests/sec you plan on having.

If the projects needs to be high performance (1000s inferences/sec) it would also be worth it to use a different (faster) lanaguage for the webserver such as Go, but that would be more complex.
