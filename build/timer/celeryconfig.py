from kombu import Exchange, Queue

BROKER_URL = 'amqp://guest@localhost//'
CELERY_RESULT_BACKEND = 'amqp://'

CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
CELERY_ACCEPT_CONTENT=['json']
CELERY_TIMEZONE = 'US/Pacific'
CELERY_ENABLE_UTC = True

CELERYD_CONCURRENCY = 1

CELERYD_PREFETCH_MULTIPLIER = 1

CELERY_DEFAULT_QUEUE = 'pisces_timer'
CELERY_QUEUES = (
    Queue ('pisces_timer', Exchange ('pisces_timer'), routing_key = 'pisces_timer'),
)