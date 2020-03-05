from fastapi import FastAPI
from starlette.responses import PlainTextResponse
from starlette.staticfiles import StaticFiles
from typing import Dict, Optional as Opt
from pydantic import BaseModel
from starlette.requests import Request

"""
* Have programme create the srimdata and pickle/sqlite it
* create api points with fast api
* uppcikle data, serve to api
* display in webpage starlette: list page and item pages
"""
"""
* buttons / inputs to run new srims
* login auth
"""


class Item(BaseModel):
    name: str
    description: Opt[str] = None
    price: float
    tax: Opt[float] = None


app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")


@app.post("/items/")
async def create_item(item: Item) -> Item:
    return item


@app.get("/items/{item_id}")
async def read_item(item_id: int) -> Dict[str, int]:
    return {"item_id": item_id}


@app.get("/")
async def home(request: Request) -> PlainTextResponse:
    return PlainTextResponse(f"Hello, world!")
